#!/usr/bin/env python3


import argparse
from typing import Any, Dict

import agglovar
import polars as pl
import pysam
from agglovar.pairwise.overlap import PairwiseOverlap


def get_chrom_sort_key(chrom: str) -> tuple:
    chrom_str = str(chrom).lower()

    if chrom_str.startswith("chr"):
        chrom_part = chrom_str[3:]
    else:
        chrom_part = chrom_str

    if chrom_part.isdigit():
        return (0, int(chrom_part), "")
    elif chrom_part == "x":
        return (1, 0, "")
    elif chrom_part == "y":
        return (1, 1, "")
    elif chrom_part in ["m", "mt"]:
        return (1, 2, "")
    else:
        return (2, 0, chrom_part)


def vcf_to_df(vcf_path: str) -> pl.DataFrame:
    variants = []
    with pysam.VariantFile(vcf_path) as vcf_in:
        header_info_keys = set(vcf_in.header.info.keys())

        for record in vcf_in:
            pos_0based = record.pos - 1

            svlen_val = (
                record.info.get("SVLEN") if "SVLEN" in header_info_keys else None
            )
            if isinstance(svlen_val, (list, tuple)):
                varlen = svlen_val[0] if svlen_val else None
            else:
                varlen = svlen_val

            end_val = record.stop
            if end_val:
                end_0based = end_val - 1
            elif varlen is not None:
                end_0based = pos_0based + abs(varlen)
            elif record.alts and not record.alts[0].startswith("<"):
                varlen = len(record.alts[0]) - len(record.ref)
                end_0based = pos_0based + abs(varlen)
            else:
                end_0based = pos_0based + 1

            if varlen is None and record.alts and not record.alts[0].startswith("<"):
                varlen = len(record.alts[0]) - len(record.ref)

            svtype = None
            if "SVTYPE" in header_info_keys:
                svtype = record.info.get("SVTYPE")
            if not svtype and "VARTYPE" in header_info_keys:
                svtype = record.info.get("VARTYPE")
            if not svtype:
                svtype = "SNV"

            varsubtype = (
                record.info.get("VARSUBTYPE")
                if "VARSUBTYPE" in header_info_keys
                else None
            )
            var_score = (
                record.info.get("SCORE") if "SCORE" in header_info_keys else None
            )

            var_data: Dict[str, Any] = {
                "chrom": record.chrom,
                "pos": pos_0based,
                "end": end_0based,
                "id": record.id if record.id else ".",
                "vartype": svtype,
                "varsubtype": varsubtype,
                "varlen": int(varlen) if varlen is not None else None,
                "var_score": float(var_score) if var_score is not None else None,
                "ref": record.ref if not record.ref.startswith("<") else "N",
                "alt": record.alts[0] if record.alts else svtype,
                "seq": (
                    record.alts[0]
                    if record.alts and not record.alts[0].startswith("<")
                    else None
                ),
            }
            variants.append(var_data)

    if not variants:
        return pl.DataFrame(schema=agglovar.schema.VARIANT)

    df = pl.DataFrame(variants, infer_schema_length=None)
    for col, dtype in agglovar.schema.VARIANT.items():
        if col not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=dtype).alias(col))
        else:
            df = df.with_columns(pl.col(col).cast(dtype, strict=False))

    return df.select(list(agglovar.schema.VARIANT.keys()))


def df_to_vcf(df: pl.DataFrame, template_vcf_path: str, output_vcf_path: str) -> None:
    df_sorted = df.to_dicts()
    df_sorted.sort(key=lambda x: (get_chrom_sort_key(x["chrom"]), x["pos"], x["end"]))

    with pysam.VariantFile(template_vcf_path) as vcf_in:
        header = vcf_in.header.copy()

        if "VARTYPE" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type">'
            )
        if "VARSUBTYPE" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=VARSUBTYPE,Number=1,Type=String,Description="Variant subtype">'
            )
        if "SVTYPE" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
            )
        if "SVLEN" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
            )
        if "END" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
            )
        if "SCORE" not in [h.name for h in header.info.values()]:
            header.add_line(
                '##INFO=<ID=SCORE,Number=1,Type=Float,Description="Variant score">'
            )

        for chrom in df.get_column("chrom").unique().to_list():
            if chrom not in header.contigs:
                header.add_line(f"##contig=<ID={chrom},length=249250621>")

    with pysam.VariantFile(output_vcf_path, "w", header=header) as vcf_out:
        for row in df_sorted:
            pos_1based = row["pos"] + 1
            end_1based = row["end"] + 1

            ref_allele = row.get("ref", "N")
            alt_allele = row.get("alt")
            vartype = row.get("vartype", "SNV")

            if not alt_allele or (
                isinstance(alt_allele, str) and alt_allele == vartype
            ):
                alt_allele = f"<{vartype}>"

            if ref_allele.startswith("<") or not ref_allele:
                ref_allele = "N"

            record = header.new_record(
                contig=row["chrom"],
                start=pos_1based - 1,
                stop=end_1based,
                alleles=(ref_allele, alt_allele),
                id=row.get("id", "."),
                qual=None,
                filter="PASS",
            )

            record.info["VARTYPE"] = vartype
            if vartype not in ["SNV", "DEL", "INS"]:
                record.info["SVTYPE"] = vartype

            varsubtype = row.get("varsubtype")
            if varsubtype:
                record.info["VARSUBTYPE"] = varsubtype

            varlen = row.get("varlen")
            if varlen is not None and abs(varlen) > 0:
                record.info["SVLEN"] = abs(varlen)

            record.stop = end_1based

            var_score = row.get("var_score")
            if var_score is not None:
                record.info["SCORE"] = float(var_score)

            vcf_out.write(record)


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge multiple VCFs using Agglovar")
    parser.add_argument("--vcfs", nargs="+", required=True, help="Input VCF files")
    parser.add_argument("--out_vcf", required=True, help="Output VCF file path")

    parser.add_argument("--ro_min", type=float, default=None)
    parser.add_argument("--size_ro_min", type=float, default=None)
    parser.add_argument("--offset_max", type=int, default=None)
    parser.add_argument("--offset_prop_max", type=float, default=None)
    parser.add_argument("--match_ref", action="store_true", default=False)
    parser.add_argument("--match_alt", action="store_true", default=False)
    parser.add_argument("--match_prop_min", type=float, default=None)

    args = parser.parse_args()

    df_cumulative = vcf_to_df(args.vcfs[0])
    join_strategy = PairwiseOverlap(
        ro_min=args.ro_min,
        size_ro_min=args.size_ro_min,
        offset_max=args.offset_max,
        offset_prop_max=args.offset_prop_max,
        match_ref=args.match_ref,
        match_alt=args.match_alt,
        match_prop_min=args.match_prop_min,
    )

    for idx, vcf_path in enumerate(args.vcfs[1:], 1):
        df_next = vcf_to_df(vcf_path)
        join_table = join_strategy.join(df_cumulative.lazy(), df_next.lazy()).collect()
        df_next_with_index = df_next.with_row_index("index")

        unmerged_next = df_next_with_index.join(
            join_table.select("index_b"),
            left_on="index",
            right_on="index_b",
            how="anti",
        ).drop("index")

        df_cumulative = pl.concat([df_cumulative, unmerged_next])

    df_to_vcf(df_cumulative, args.vcfs[0], args.out_vcf)


if __name__ == "__main__":
    main()
