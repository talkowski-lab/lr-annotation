#!/usr/bin/env python3


import argparse
from pathlib import Path
from typing import Optional
from datetime import datetime

import polars as pl
import pysam


def process_parquet_to_vcf(
    parquet_file: Path,
    ref_fa: Path,
    limit: Optional[int]
) -> pl.DataFrame:
    df = pl.read_parquet(parquet_file)
    if limit:
        df = df.head(limit)
    fasta = pysam.FastaFile(str(ref_fa))

    def parse_ref(row):
        vartype = row.get('vartype', 'SNV')

        if vartype == 'SNV' and 'ref' in df.columns and row.get('ref') is not None:
            return row['ref'].upper()
        elif vartype == 'DEL':
            ref_base = fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1)
            if 'seq' in df.columns and row.get('seq') is not None:
                deleted_seq = row['seq']
            else:
                deleted_seq = fasta.fetch(row['chrom'], row['pos'] + 1, row['end'])
            return (ref_base + deleted_seq).upper()
        elif vartype == 'INS':
            return fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1).upper()
        else:
            return fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1).upper()

    def parse_alt(row):
        vartype = row.get('vartype', 'SNV')

        if vartype == 'SNV' and 'alt' in df.columns and row.get('alt') is not None:
            return row['alt'].upper()
        elif vartype == 'DEL':
            return fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1).upper()
        elif vartype == 'INS':
            ref_base = fasta.fetch(row['chrom'], row['pos'], row['pos'] + 1)
            if 'seq' in df.columns and row.get('seq') is not None:
                inserted_seq = row['seq']
            else:
                inserted_seq = ''
            return (ref_base + inserted_seq).upper()
        else:
            return f"<{vartype}>".upper()

    def parse_info(row):
        info_fields = []
        
        vartype = row.get('vartype', 'SNV')
        if vartype:
            info_fields.append(f"VARTYPE={vartype}")
            if vartype not in ['SNV', 'DEL', 'INS']:
                info_fields.append(f"SVTYPE={vartype}")
        
        varsubtype = row.get('varsubtype')
        if varsubtype:
            info_fields.append(f"VARSUBTYPE={varsubtype}")
        
        if 'varlen' in df.columns and row.get('varlen') is not None:
            svlen = abs(row['varlen'])
        elif 'pos' in df.columns and 'end' in df.columns and row.get('pos') is not None and row.get('end') is not None:
            svlen = abs(row['end'] - row['pos'])
        else:
            svlen = None
        
        if svlen is not None and svlen > 0:
            info_fields.append(f"SVLEN={svlen}")
        
        end = row.get('end')
        if end is not None:
            info_fields.append(f"END={end + 1}")
        
        var_score = row.get('var_score')
        if var_score is not None:
            info_fields.append(f"SCORE={var_score:.4f}")
        
        return ';'.join(info_fields) if info_fields else '.'

    rows = []
    for row_dict in df.iter_rows(named=True):
        vcf_pos = row_dict['pos'] + 1
        row_values = [
            row_dict['chrom'],
            vcf_pos,
            row_dict['id'],
            parse_ref(row_dict),
            parse_alt(row_dict),
            '.',
            'PASS',
            parse_info(row_dict),
            'GT',
            '0/1'
        ]
        rows.append(row_values)

    fasta.close()
    schema = {
        'CHROM': pl.Utf8,
        'POS': pl.Int64,
        'ID': pl.Utf8,
        'REF': pl.Utf8,
        'ALT': pl.Utf8,
        'QUAL': pl.Utf8,
        'FILTER': pl.Utf8,
        'INFO': pl.Utf8,
        'FORMAT': pl.Utf8,
        'SAMPLE': pl.Utf8
    }
    vcf_df = pl.DataFrame(rows, schema=schema, orient='row')
    return vcf_df


def write_vcf_header(f, sample_id: str, ref_fa: Path, parquet_files: list):
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
    f.write(f"##reference={ref_fa}\n")
    f.write(f"##source=PAV3_parquet_to_vcf\n")
    
    fasta = pysam.FastaFile(str(ref_fa))
    for contig in fasta.references:
        length = fasta.get_reference_length(contig)
        f.write(f"##contig=<ID={contig},length={length}>\n")
    fasta.close()
    
    f.write('##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type">\n')
    f.write('##INFO=<ID=VARSUBTYPE,Number=1,Type=String,Description="Variant subtype">\n')
    f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    f.write('##INFO=<ID=SCORE,Number=1,Type=Float,Description="Variant score">\n')
    f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')


def main():
    parser = argparse.ArgumentParser(
        description="Convert PAV parquet files to VCF format",
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
        help='Output VCF file path'
    )
    parser.add_argument(
        '--ref',
        required=True,
        help='Reference FASTA file (for extracting ref/alt bases)'
    )
    parser.add_argument(
        '--sample',
        required=True,
        help='Sample ID for VCF'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Output only the first N rows (optional)'
    )

    args = parser.parse_args()

    parquet_files = [Path(p) for p in args.parquet]
    output_vcf = Path(args.output)
    ref_fa = Path(args.ref)
    sample_id = args.sample

    dfs = []
    for parquet_file in parquet_files:
        df = process_parquet_to_vcf(parquet_file, ref_fa, args.limit)
        dfs.append(df)
    combined_df = pl.concat(dfs)

    combined_df = (combined_df
        .sort('SAMPLE', descending=True, nulls_last=True)
        .unique(subset=['CHROM', 'POS', 'REF', 'ALT'], keep='first')
    )

    chr_order = {c: i for i, c in enumerate([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])}
    combined_df = (combined_df
        .with_columns(pl.col('CHROM').replace_strict(chr_order, default=999).alias('_chr_order'))
        .sort(['_chr_order', 'CHROM', 'POS'])
        .drop('_chr_order')
    )

    combined_df = combined_df.rename({'SAMPLE': sample_id})

    with open(output_vcf, 'w') as f:
        write_vcf_header(f, sample_id, ref_fa, parquet_files)
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")
        combined_df.write_csv(f, separator='\t', include_header=False, null_value='.')


if __name__ == '__main__':
    main()
