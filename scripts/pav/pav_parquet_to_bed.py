#!/usr/bin/env python3


import argparse
import sys
from pathlib import Path
from typing import List, Optional, Set

import polars as pl
import pysam


def get_all_columns(parquet_files: List[Path]) -> List[str]:
    all_columns = set()
    for pf in parquet_files:
        df = pl.read_parquet(pf)
        all_columns.update(df.columns)
    
    base_columns = ['chrom', 'pos', 'end', 'id', 'ref', 'alt']
    other_columns = sorted([col for col in all_columns if col not in base_columns])
    
    return base_columns + other_columns


def process_parquet_to_bed(
    parquet_file: Path, 
    all_columns: List[str],
    ref_fa: Optional[Path] = None,
    limit: Optional[int] = None
) -> pl.DataFrame:
    df = pl.read_parquet(parquet_file)
    if limit:
        df = df.head(limit)
    fasta = pysam.FastaFile(str(ref_fa))
    
    def get_ref_column(row):
        if 'ref' in df.columns and row.get('ref') is not None:
            return str(row['ref']).upper()
        elif 'seq' in df.columns and row.get('seq') is not None:
            if row.get('vartype') == 'DEL':
                return fasta.fetch(row['chrom'], row['pos'] - 1, row['end']).upper()
            else:
                return fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1).upper()
        else:
            return 'N'
    
    def get_alt_column(row):
        if 'alt' in df.columns and row.get('alt') is not None:
            return str(row['alt']).upper()
        elif 'seq' in df.columns and row.get('seq') is not None:
            if row.get('vartype') == 'DEL':
                return fasta.fetch(row['chrom'], row['pos'] - 1, row['pos']).upper()
            else:
                ref_base = fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1)
                return (ref_base + row['seq']).upper()
        elif 'vartype' in df.columns and row.get('vartype') is not None:
            return f"<{row['vartype']}>"
        else:
            return '<SNV>'
    
    rows = []
    for row_dict in df.iter_rows(named=True):
        ref_val = get_ref_column(row_dict)
        alt_val = get_alt_column(row_dict)
        
        variant_id = row_dict.get('id', 'Missing_ID') if row_dict.get('id') is not None else 'Missing_ID'
        
        row_values = [
            str(row_dict['chrom']),
            row_dict['pos'],
            row_dict['end'],
            str(variant_id),
            str(ref_val),
            str(alt_val)
        ]
        
        for col in all_columns[6:]:
            val = row_dict.get(col)
            if val is None:
                row_values.append(None)
            elif isinstance(val, (list, dict)):
                row_values.append(str(val))
            else:
                row_values.append(str(val))
        
        rows.append(row_values)
    
    if fasta:
        fasta.close()
    
    schema = {
        'chrom': pl.Utf8,
        'pos': pl.Int64,
        'end': pl.Int64,
        'id': pl.Utf8,
        'ref': pl.Utf8,
        'alt': pl.Utf8
    }
    for col in all_columns[6:]:
        schema[col] = pl.Utf8
    
    bed_df = pl.DataFrame(rows, schema=schema, orient='row')
    
    return bed_df


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
    parser.add_argument(
        '--ref',
        required=True,
        help='Reference FASTA file (for extracting ref/alt bases)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Output only the first N rows (optional)'
    )
    
    args = parser.parse_args()
    
    parquet_files = [Path(p) for p in args.parquet]
    output_bed = Path(args.output)
    ref_fa = Path(args.ref) if args.ref else None
    
    all_columns = get_all_columns(parquet_files)
    
    dfs = []
    for parquet_file in parquet_files:
        df = process_parquet_to_bed(parquet_file, all_columns, ref_fa, args.limit)
        dfs.append(df)
    
    combined_df = pl.concat(dfs)
    
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    combined_df = combined_df.with_columns(
        pl.when(pl.col('chrom').is_in(chr_order))
        .then(pl.col('chrom').replace_strict(
            {c: i for i, c in enumerate(chr_order)},
            default=999
        ))
        .otherwise(999)
        .alias('_chr_order')
    )
    
    combined_df = combined_df.sort(['_chr_order', 'chrom', 'pos', 'end']).drop('_chr_order')
    
    with open(output_bed, 'w') as f:
        combined_df.write_csv(f, separator='\t', include_header=True, null_value='')

if __name__ == '__main__':
    main()
