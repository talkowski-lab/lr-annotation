#!/usr/bin/env python3

import argparse
import logging
import sys
from typing import Any, Dict

import agglovar
import polars as pl
import pysam
from agglovar.join.pair import PairwiseIntersect

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)


def vcf_to_df(vcf_path: str) -> pl.DataFrame:
    logging.info(f"Converting VCF to DataFrame: {vcf_path}")
    
    variants = []
    with pysam.VariantFile(vcf_path) as vcf_in:
        for record in vcf_in:
            varlen = record.info.get("SVLEN", [0])[0]
            if not varlen and record.alts:
                varlen = len(record.alts[0]) - len(record.ref)

            end = record.info.get("END", record.pos + abs(varlen))
            if record.info.get("SVTYPE") == "INS" and not record.info.get("END"):
                end = record.pos + 1

            var_data: Dict[str, Any] = {
                "chrom": record.chrom,
                "pos": record.pos,
                "end": end,
                "id": record.id,
                "vartype": record.info.get("SVTYPE", "UNKNOWN"),
                "varlen": varlen,
                "ref": record.ref,
                "alt": record.alts[0] if record.alts else None,
                "seq": (
                    record.alts[0] if record.alts and ">" not in record.alts[0] else None
                ),
            }
            variants.append(var_data)

    df = pl.DataFrame(variants)
    logging.info(f"Parsed {len(variants)} variants from {vcf_path}")
    
    for col, dtype in agglovar.schema.VARIANT.items():
        if col not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=dtype).alias(col))
        else:
            df = df.with_columns(pl.col(col).cast(dtype, strict=False))

    logging.info(f"Schema conversion complete for {vcf_path}")
    return df.select(list(agglovar.schema.VARIANT.keys()))


def df_to_vcf(df: pl.DataFrame, template_vcf_path: str, output_vcf_path: str) -> None:
    logging.info(f"Writing {len(df)} variants to {output_vcf_path}")
    with pysam.VariantFile(template_vcf_path) as vcf_in:
        header = vcf_in.header
        for chrom in df.get_column("chrom").unique().to_list():
            if chrom not in header.contigs:
                header.add_line(f"##contig=<ID={chrom},length=249250621>")

    with pysam.VariantFile(output_vcf_path, "w", header=header) as vcf_out:
        for row in df.iter_rows(named=True):
            pos = int(row["pos"])
            end = int(row["end"])

            record = header.new_record(
                contig=row["chrom"],
                start=pos - 1,
                stop=end,
                alleles=(row["ref"], row["alt"]),
                id=row["id"],
            )

            record.info.clear()
            if row["vartype"]:
                record.info["SVTYPE"] = row["vartype"]
            if row["varlen"]:
                record.info["SVLEN"] = int(row["varlen"])
            if row["vartype"] != "INS":
                record.info["END"] = end

            vcf_out.write(record)


def main() -> None:
    logging.info("Starting AggloVar merge script")
    
    parser = argparse.ArgumentParser(
        description="Merge multiple VCFs using AggloVar"
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
    
    logging.info(f"Input VCFs: {args.vcfs}")
    logging.info(f"Output VCF: {args.out_vcf}")

    if not args.vcfs:
        raise ValueError("At least one input VCF is required.")

    logging.info(f"Loading first VCF: {args.vcfs[0]}")
    df_cumulative = vcf_to_df(args.vcfs[0])
    logging.info(f"Loaded {len(df_cumulative)} variants from {args.vcfs[0]}")

    logging.info("Initializing PairwiseOverlap strategy")
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
        logging.info(f"Processing VCF {idx}/{len(args.vcfs)-1}: {vcf_path}")
        df_next = vcf_to_df(vcf_path)
        logging.info(f"Loaded {len(df_next)} variants from {vcf_path}")

        if len(df_next) == 0:
            logging.warning(f"Skipping empty VCF: {vcf_path}")
            continue

        logging.info(f"Running pairwise merge with {len(df_cumulative)} cumulative variants")
        join_table = join_strategy.join(df_cumulative.lazy(), df_next.lazy()).collect()
        logging.info(f"Found {len(join_table)} overlapping variants")

        df_next_with_index = df_next.with_row_index("index")

        unmerged_next = df_next_with_index.join(
            join_table.select("index_b"),
            left_on="index",
            right_on="index_b",
            how="anti",
        ).drop("index")

        logging.info(f"Adding {len(unmerged_next)} new variants to cumulative set")
        df_cumulative = pl.concat([df_cumulative, unmerged_next])
        logging.info(f"Cumulative set now has {len(df_cumulative)} variants")

    logging.info(f"Total variants in final merged set: {len(df_cumulative)}")
    df_to_vcf(df_cumulative, args.vcfs[0], args.out_vcf)
    logging.info("Merge complete!")


if __name__ == "__main__":
  main()
