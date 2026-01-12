#!/usr/bin/env python3


import argparse
from pathlib import Path

import pysam


def combine_haplotypes(input_vcf: Path, output_vcf: Path, sample_name: str) -> None:
    vcf_in = pysam.VariantFile(str(input_vcf))
    input_samples = list(vcf_in.header.samples)
    mat_sample = input_samples[0]
    pat_sample = input_samples[1]

    header_out = pysam.VariantHeader()
    for record in vcf_in.header.records:
        if record.key != "SAMPLE":
            header_out.add_record(record)
    header_out.add_sample(sample_name)
    vcf_out = pysam.VariantFile(str(output_vcf), "w", header=header_out)

    for rec in vcf_in:
        rec_out = vcf_out.new_record(
            contig=rec.chrom,
            start=rec.pos - 1,
            stop=rec.stop,
            alleles=tuple([rec.ref] + list(rec.alts)) if rec.alts else (rec.ref,),
            id=rec.id,
            qual=rec.qual,
            filter=rec.filter,
        )

        for key, val in rec.info.items():
            rec_out.info[key] = val

        gt_mat = rec.samples[mat_sample]["GT"]
        gt_pat = rec.samples[pat_sample]["GT"]

        mat_has_alt = any(g is not None and g > 0 for g in gt_mat)
        pat_has_alt = any(g is not None and g > 0 for g in gt_pat)

        if mat_has_alt and pat_has_alt:
            rec_out.samples[sample_name]["GT"] = (1, 1)
        elif mat_has_alt or pat_has_alt:
            rec_out.samples[sample_name]["GT"] = (0, 1)
        else:
            rec_out.samples[sample_name]["GT"] = (0, 0)

        vcf_out.write(rec_out)

    vcf_in.close()
    vcf_out.close()


def main():
    parser = argparse.ArgumentParser(
        description="Combine maternal and paternal haplotype VCFs into a diploid VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file with two haplotype samples (maternal and paternal)",
    )
    parser.add_argument("--output", required=True, help="Output diploid VCF file")
    parser.add_argument(
        "--sample", required=True, help="Sample name for output diploid VCF"
    )

    args = parser.parse_args()

    input_vcf = Path(args.input)
    output_vcf = Path(args.output)

    combine_haplotypes(input_vcf, output_vcf, args.sample)


if __name__ == "__main__":
    main()
