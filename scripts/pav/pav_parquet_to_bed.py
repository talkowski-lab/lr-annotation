#!/usr/bin/env python3


import argparse
import sys
from pathlib import Path
from typing import List

import polars as pl


def process_parquet_to_bed(parquet_file: Path, output_bed: Path, append: bool = False) -> None:
    df = pl.read_parquet(parquet_file)
    
    columns_to_select = [
        pl.col('chrom'),
        pl.col('pos'),
        pl.col('end'),
        pl.col('id'),
    ]
    
    if 'vartype' in df.columns:
        columns_to_select.append(pl.col('vartype').fill_null('SNV'))
    else:
        columns_to_select.append(pl.lit('SNV').alias('vartype'))
    
    bed_df = df.select(columns_to_select)
    mode = 'a' if append else 'w'
    with open(output_bed, mode) as f:
        bed_df.write_csv(f, separator='\t', include_header=False)


def main():
    parser = argparse.ArgumentParser(
        description="Convert PAV parquet files to BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--parquet',
        nargs='+',
        required=True,
        help='PAV parquet file(s) to convert'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output BED file path'
    )
    
    args = parser.parse_args()
    
    parquet_files = [Path(p) for p in args.parquet]
    output_bed = Path(args.output)
    
    for i, parquet_file in enumerate(parquet_files):
        append = i > 0
        process_parquet_to_bed(parquet_file, output_bed, append)

if __name__ == '__main__':
    main()
