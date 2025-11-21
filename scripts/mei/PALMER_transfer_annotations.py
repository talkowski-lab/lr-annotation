#!/usr/bin/env python3


import argparse
import gzip
import sys
from difflib import SequenceMatcher

from Bio import SeqIO
from pysam import VariantFile


def parse_intersection_line(line: str):
    fields = line.strip().split("\t")
    return {
        "chrom": fields[0],
        "sv_positions": fields[3].split(","),
        "sv_refs": fields[4].split(","),
        "sv_lengths": fields[5].split(","),
        "mei_start_0based": int(fields[7]),
        "mei_pos_1based": int(fields[8]),
        "mei_id": fields[9],
        "mei_length": abs(int(fields[10])),
        "mei_samples": fields[11].rstrip(",").split(","),
    }


def compute_size_similarity(target_length: int, mei_length: int) -> float:
    if target_length <= 0 or mei_length <= 0:
        return 0.0
    return min(target_length, mei_length) / max(target_length, mei_length)


def compute_reciprocal_overlap(
    target_start: int, target_length: int, mei_start: int, mei_length: int
) -> float:
    target_end = target_start + target_length
    mei_end = mei_start + mei_length
    overlap = min(target_end, mei_end) - max(target_start, mei_start)
    if overlap <= 0:
        return 0.0
    return min(overlap / target_length, overlap / mei_length)


def compute_sequence_similarity(seq_a: str, seq_b: str) -> float:
    if not seq_a and not seq_b:
        return 1.0
    if not seq_a or not seq_b:
        return 0.0
    if seq_a == seq_b:
        return 1.0
    return SequenceMatcher(None, seq_a, seq_b).ratio()


def record_passes_filters(
    record,
    ref_allele: str,
    sv_length: int,
    mei_data: dict,
    inserted_sequence: str,
    samples: list[str],
    args: argparse.Namespace,
) -> bool:
    svlen_field = record.info.get("SVLEN")
    svlen = abs(svlen_field[0]) if isinstance(svlen_field, tuple) else abs(svlen_field)
    target_length = abs(sv_length)
    if svlen != target_length - 1:
        return False

    svtype_field = record.info.get("SVTYPE")
    svtype = svtype_field[0] if isinstance(svtype_field, tuple) else svtype_field
    if svtype != "INS":
        return False

    if record.ref != ref_allele:
        return False

    if not record.alts:
        return False

    alt_allele = record.alts[0]
    sequence_similarity = compute_sequence_similarity(alt_allele, inserted_sequence)
    if sequence_similarity < args.sequence_similarity:
        return False

    size_similarity = compute_size_similarity(target_length, mei_data["mei_length"])
    if size_similarity < args.size_similarity:
        return False

    reciprocal_overlap = compute_reciprocal_overlap(
        record.pos - 1, target_length, mei_data["mei_start_0based"], mei_data["mei_length"]
    )
    if reciprocal_overlap < args.reciprocal_overlap:
        return False

    breakpoint_distance = abs(record.pos - mei_data["mei_pos_1based"])
    if breakpoint_distance > args.breakpoint_window:
        return False

    genotypes = [record.samples[sample]["GT"] for sample in samples]
    sv_samples = [sample for gt, sample in zip(genotypes, samples) if 1 in gt]
    shared = set(sv_samples) & set(mei_data["mei_samples"])
    if len(shared) < args.min_shared_samples:
        return False

    return True


def write_header(handle):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write(
        "##INFO=<ID=ME_TYPE,Number=1,Type=String,Description=\"Type of mobile element\">\n"
    )
    handle.write(
        "##INFO=<ID=ME_LEN,Number=1,Type=Integer,Description=\"Length of PALMER call\">\n"
    )
    handle.write(
        "##INFO=<ID=ME_ID,Number=1,Type=String,Description=\"ID of the PALMER variant\">\n"
    )
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Transfer PALMER MEI annotations onto a short-read VCF"
    )
    parser.add_argument("--intersection", required=True, help="Path to bedtools intersect output")
    parser.add_argument("--target-vcf", required=True, help="Target VCF to annotate")
    parser.add_argument("--ins-fa", required=True, help="FASTA with MEI insert sequences")
    parser.add_argument("--me-type", required=True, help="MEI class under evaluation")
    parser.add_argument(
        "--output",
        required=True,
        help="Destination bgzipped VCF containing transferable annotations",
    )
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        required=True,
        help="Minimum reciprocal overlap between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--size-similarity",
        type=float,
        required=True,
        help="Minimum size similarity between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--sequence-similarity",
        type=float,
        required=True,
        help="Minimum sequence similarity between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--breakpoint-window",
        type=int,
        required=True,
        help="Maximum distance (bp) between event breakpoints",
    )
    parser.add_argument(
        "--min-shared-samples",
        type=int,
        required=True,
        help="Minimum number of shared samples required for annotation",
    )
    args = parser.parse_args()

    target_vcf = VariantFile(args.target_vcf)
    ins_fasta = SeqIO.index(args.ins_fa, "fasta")
    samples = list(target_vcf.header.samples)

    with open(args.intersection, "r", encoding="utf-8") as source, gzip.open(args.output, "wt", encoding="utf-8") as output_handle:
        write_header(output_handle)
        for raw_line in source:
            mei_data = parse_intersection_line(raw_line)
            chrom = mei_data["chrom"]
            for ref, pos, length in zip(mei_data["sv_refs"], mei_data["sv_positions"], mei_data["sv_lengths"]):
                sv_pos = int(pos)
                sv_length = int(length)
                ref_allele = ref.split("_")[0]
                sv_id = f"{chrom}:{pos};{ref}"
                inserted_sequence = str(ins_fasta[sv_id].seq)

                for record in target_vcf.fetch(chrom, sv_pos - 1, sv_pos):
                    if "ME_TYPE" in record.info:
                        continue

                    if not record_passes_filters(record, ref_allele, sv_length, mei_data,
                                                 inserted_sequence, samples, args):
                        continue

                    output_handle.write(
                        f"{chrom}\t{sv_pos}\t.\t{record.ref}\t{record.alts[0]}\t.\t.\t"
                    )
                    output_handle.write(
                        f"ME_TYPE={args.me_type};ME_LEN={mei_data['mei_length']};ME_ID={mei_data['mei_id']}\n"
                    )

    target_vcf.close()
    if hasattr(ins_fasta, "close"):
        ins_fasta.close()


if __name__ == "__main__":
    main()
