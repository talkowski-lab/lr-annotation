#!/usr/bin/env python3


import argparse
from pathlib import Path
from typing import Optional

import polars as pl
import pysam


def process_parquet_to_bed(
    parquet_file: Path, ref_fa: Path, limit: Optional[int]
) -> pl.DataFrame:
    df = pl.read_parquet(parquet_file)
    if limit:
        df = df.head(limit)
    fasta = pysam.FastaFile(str(ref_fa))

    def parse_ref(row):
        vartype = row.get("vartype", "SNV")

        if vartype == "SNV" and "ref" in df.columns and row.get("ref") is not None:
            return row["ref"].upper()
        elif vartype == "DEL":
            ref_base = fasta.fetch(row["chrom"], row["pos"], row["pos"] + 1)
            if "seq" in df.columns and row.get("seq") is not None:
                deleted_seq = row["seq"]
            else:
                deleted_seq = fasta.fetch(row["chrom"], row["pos"] + 1, row["end"])
            return (ref_base + deleted_seq).upper()
        elif vartype == "INS":
            return fasta.fetch(row["chrom"], row["pos"], row["pos"] + 1).upper()
        else:
            return "N"

    def parse_alt(row):
        vartype = row.get("vartype", "SNV")

        if vartype == "SNV" and "alt" in df.columns and row.get("alt") is not None:
            return row["alt"].upper()
        elif vartype == "DEL":
            return fasta.fetch(row["chrom"], row["pos"], row["pos"] + 1).upper()
        elif vartype == "INS":
            ref_base = fasta.fetch(row["chrom"], row["pos"], row["pos"] + 1)
            if "seq" in df.columns and row.get("seq") is not None:
                inserted_seq = row["seq"]
            else:
                inserted_seq = ""
            return (ref_base + inserted_seq).upper()
        else:
            return f"<{vartype}>".upper()

    def parse_svlen(row):
        if "varlen" in df.columns and row.get("varlen") is not None:
            return abs(row["varlen"])
        elif (
            "pos" in df.columns
            and "end" in df.columns
            and row.get("pos") is not None
            and row.get("end") is not None
        ):
            return abs(row["end"] - row["pos"])
        else:
            return -1

    rows = []
    for row_dict in df.iter_rows(named=True):
        row_values = [
            row_dict["chrom"],
            row_dict["pos"],
            row_dict["end"],
            row_dict["id"],
            parse_ref(row_dict),
            parse_alt(row_dict),
            parse_svlen(row_dict),
            row_dict.get("vartype", "SNV"),
            row_dict.get("varsubtype", None),
            row_dict.get("var_score", None),
        ]
        rows.append(row_values)

    fasta.close()
    schema = {
        "chrom": pl.Utf8,
        "pos": pl.Int64,
        "end": pl.Int64,
        "id": pl.Utf8,
        "ref": pl.Utf8,
        "alt": pl.Utf8,
        "svlen": pl.Int64,
        "vartype": pl.Utf8,
        "varsubtype": pl.Utf8,
        "score": pl.Float64,
    }
    bed_df = pl.DataFrame(rows, schema=schema, orient="row")
    return bed_df


def main():
    parser = argparse.ArgumentParser(
        description="Convert PAV parquet files to BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--parquet", nargs="+", required=True, help="PAV parquet file(s) to convert"
    )
    parser.add_argument("--output", required=True, help="Output BED file path")
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file (for extracting ref/alt bases)",
    )
    parser.add_argument(
        "--limit", type=int, help="Output only the first N rows (optional)"
    )

    args = parser.parse_args()

    parquet_files = [Path(p) for p in args.parquet]
    output_bed = Path(args.output)
    ref_fa = Path(args.ref) if args.ref else None

    dfs = []
    for parquet_file in parquet_files:
        df = process_parquet_to_bed(parquet_file, ref_fa, args.limit)
        dfs.append(df)
    combined_df = pl.concat(dfs)

    combined_df = combined_df.sort("score", descending=True, nulls_last=True).unique(
        subset=["chrom", "pos", "end", "ref", "alt"], keep="first"
    )

    chr_order = {
        c: i
        for i, c in enumerate(
            [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        )
    }
    combined_df = (
        combined_df.with_columns(
            pl.col("chrom").replace_strict(chr_order, default=999).alias("_chr_order")
        )
        .sort(["_chr_order", "chrom", "pos", "end"])
        .drop("_chr_order")
    )

    with open(output_bed, "w") as f:
        combined_df.write_csv(f, separator="\t", include_header=True, null_value="")


if __name__ == "__main__":
    main()
