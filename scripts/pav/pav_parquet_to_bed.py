#!/usr/bin/env python3


import argparse
import sys
from pathlib import Path
from typing import List, Optional

import polars as pl
import pysam


def process_parquet_to_bed(
    parquet_file: Path, 
    output_bed: Path, 
    ref_fa: Optional[Path] = None,
    append: bool = False
) -> None:
    df = pl.read_parquet(parquet_file)
    fasta = pysam.FastaFile(str(ref_fa))
    
    def get_ref_base(row):
        try:
            return fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1)
        except:
            return '.'
    
    def get_ref_column(row):
        if 'ref' in df.columns and row.get('ref') is not None:
            return str(row['ref']).upper()
        elif 'seq' in df.columns and row.get('seq') is not None:
            return get_ref_base(row).upper()
        else:
            return '.'
    
    def get_alt_column(row):
        if 'alt' in df.columns and row.get('alt') is not None:
            return str(row['alt']).upper()
        elif 'seq' in df.columns and row.get('seq') is not None:
            ref_base = get_ref_base(row)
            return (ref_base + row['seq']).upper()
        elif 'vartype' in df.columns and row.get('vartype') is not None:
            return f"<{row['vartype']}>"
        else:
            return '<SNV>'
    
    rows = []
    for row_dict in df.iter_rows(named=True):
        ref_val = get_ref_column(row_dict)
        alt_val = get_alt_column(row_dict)
        
        variant_id = row_dict.get('id', '.')
        vartype_val = row_dict.get('vartype', 'SNV') if row_dict.get('vartype') is not None else 'SNV'
        
        rows.append([
            row_dict['chrom'],
            row_dict['pos'],
            row_dict['end'],
            ref_val,
            alt_val,
            variant_id,
            vartype_val
        ])
    
    if fasta:
        fasta.close()
    
    bed_df = pl.DataFrame(rows, schema=['chrom', 'pos', 'end', 'ref', 'alt', 'id', 'vartype'])
    
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
    parser.add_argument(
        '--ref',
        required=True,
        help='Reference FASTA file (for extracting ref/alt bases)'
    )
    
    args = parser.parse_args()
    
    parquet_files = [Path(p) for p in args.parquet]
    output_bed = Path(args.output)
    ref_fa = Path(args.ref) if args.ref else None
    
    for i, parquet_file in enumerate(parquet_files):
        append = i > 0
        process_parquet_to_bed(parquet_file, output_bed, ref_fa, append)

if __name__ == '__main__':
    main()
