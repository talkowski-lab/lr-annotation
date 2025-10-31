#!/usr/bin/env python3


import argparse
import sys
from typing import Any, Dict

import agglovar
import polars as pl
import pysam
from agglovar.join.pair import PairwiseIntersect


def get_chrom_sort_key(chrom: str) -> tuple:
    chrom_str = str(chrom).lower()
    
    if chrom_str.startswith('chr'):
        chrom_part = chrom_str[3:]
    else:
        chrom_part = chrom_str
    
    if chrom_part.isdigit():
        return (0, int(chrom_part), '')
    elif chrom_part == 'x':
        return (1, 0, '')
    elif chrom_part == 'y':
        return (1, 1, '')
    elif chrom_part in ['m', 'mt']:
        return (1, 2, '')
    else:
        return (2, 0, chrom_part)


def vcf_to_df(vcf_path: str) -> pl.DataFrame:    
    variants = []
    with pysam.VariantFile(vcf_path) as vcf_in:
        for record in vcf_in:
            svlen_val = record.info.get("SVLEN")
            if isinstance(svlen_val, (list, tuple)):
                varlen = svlen_val[0] if svlen_val else 0
            else:
                varlen = svlen_val
            
            if not varlen and record.alts:
                varlen = len(record.alts[0]) - len(record.ref)

            svtype = record.info.get("SVTYPE", "Unknown")

            end = record.stop if record.stop else record.pos + abs(varlen)

            var_data: Dict[str, Any] = {
                "chrom": record.chrom,
                "pos": record.pos,
                "end": end,
                "id": record.id,
                "vartype": svtype,
                "varlen": varlen,
                "ref": "N",
                "alt": svtype,
                "seq": record.alts[0] if record.alts else None
            }
            variants.append(var_data)

    df = pl.DataFrame(variants)
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
        header = vcf_in.header
        for chrom in df.get_column("chrom").unique().to_list():
            if chrom not in header.contigs:
                header.add_line(f"##contig=<ID={chrom},length=249250621>")

    with pysam.VariantFile(output_vcf_path, "w", header=header) as vcf_out:
        for row in df_sorted:
            record = header.new_record(
                contig=row["chrom"],
                start=row["pos"],
                stop=row["end"],
                alleles=("N", row["vartype"]),
                id=row["id"],
            )

            record.info.clear()
            if row["vartype"]:
                record.info["SVTYPE"] = row["vartype"]
            if row["varlen"]:
                record.info["SVLEN"] = row["varlen"]

            vcf_out.write(record)


def main() -> None:    
    parser = argparse.ArgumentParser(
        description="Merge multiple VCFs using Agglovar"
    )
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
    join_strategy = PairwiseIntersect(
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
