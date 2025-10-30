#!/usr/bin/env python3


import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import polars as pl
import pysam

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)


def load_pav_calls(call_files: List[Path]) -> pl.DataFrame:
    dfs = []
    for file_path in call_files:
        logging.info(f"Loading {file_path}")
        df = pl.read_parquet(file_path)
        logging.info(f"  Loaded {len(df)} variants")
        dfs.append(df)
    
    combined = pl.concat(dfs)
    return combined


def create_vcf_header(
    ref_fa: Optional[Path],
    sample_name: str,
    filters_from_data: Optional[List[str]] = None
) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    
    if ref_fa and ref_fa.exists():
        with pysam.FastaFile(str(ref_fa)) as fa:
            for chrom in fa.references:
                header.add_line(f'##contig=<ID={chrom},length={fa.get_reference_length(chrom)}>')
    
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.add_line('##INFO=<ID=PAV_FILTER,Number=.,Type=String,Description="PAV filter flags">')
    header.add_line('##INFO=<ID=PAV_SOURCE,Number=1,Type=String,Description="PAV call source">')
    header.add_line('##INFO=<ID=PAV_SCORE,Number=1,Type=Float,Description="PAV variant score">')
    header.add_line('##INFO=<ID=QRY_ID,Number=1,Type=String,Description="Query/haplotype ID">')
    header.add_line('##INFO=<ID=QRY_POS,Number=1,Type=Integer,Description="Position in query">')
    header.add_line('##INFO=<ID=QRY_END,Number=1,Type=Integer,Description="End position in query">')
    header.add_line('##INFO=<ID=QRY_REV,Number=0,Type=Flag,Description="Query is reverse complemented">')
    header.add_line('##INFO=<ID=INV_DUP,Number=0,Type=Flag,Description="Inverted duplication">')
    header.add_line('##INFO=<ID=DERIVED,Number=1,Type=String,Description="Derived from complex variant">')
    
    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    
    if filters_from_data:
        for filt in sorted(set(filters_from_data)):
            if filt and filt != 'PASS':
                header.add_line(f'##FILTER=<ID={filt},Description="Filter from PAV: {filt}">')
    
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    
    header.add_sample(sample_name)
    
    return header


def variant_to_vcf_record(
    variant: Dict,
    header: pysam.VariantHeader,
    sample_name: str
) -> pysam.VariantRecord:
    chrom = variant['chrom']
    pos = variant['pos']
    end = variant['end']
    var_id = variant['id']
    vartype = variant['vartype']
    varlen = variant['varlen']
    
    record = header.new_record(
        contig=chrom,
        start=pos,
        stop=end,
        id=var_id,
        alleles=('N', f'<{vartype}>'),
    )
    
    record.info['SVTYPE'] = vartype
    if varlen is not None:
        record.info['SVLEN'] = varlen
    
    if variant['filter'] and len(variant['filter']) > 0:
        record.info['PAV_FILTER'] = ','.join(variant['filter'])
    
    if variant['call_source']:
        record.info['PAV_SOURCE'] = variant['call_source']
    
    if variant['var_score'] is not None:
        record.info['PAV_SCORE'] = variant['var_score']
    
    if variant['qry_id']:
        record.info['QRY_ID'] = variant['qry_id']
    
    if variant['qry_pos'] is not None:
        record.info['QRY_POS'] = variant['qry_pos']
    
    if variant['qry_end'] is not None:
        record.info['QRY_END'] = variant['qry_end']
    
    if variant['qry_rev']:
        record.info['QRY_REV'] = True
    
    if variant['inv_dup']:
        record.info['INV_DUP'] = True
    
    if variant['derived']:
        record.info['DERIVED'] = variant['derived']
    
    if variant['filter'] and len(variant['filter']) > 0:
        for filt in variant['filter']:
            record.filter.add(filt)
    else:
        record.filter.add('PASS')
    
    record.samples[sample_name]['GT'] = (1,)
    
    return record


def write_vcf(
    df: pl.DataFrame,
    output_vcf: Path,
    sample_name: str,
    ref_fa: Optional[Path] = None
) -> None:
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    
    df = df.with_columns(
        pl.when(pl.col('chrom').is_in(chr_order))
        .then(pl.col('chrom').replace(
            {c: i for i, c in enumerate(chr_order)},
            default=999
        ))
        .otherwise(999)
        .alias('_chr_order')
    )
    
    df = df.sort(['_chr_order', 'chrom', 'pos']).drop('_chr_order')
    
    filters_from_data = df.get_column('filter').explode().unique().drop_nulls().to_list()
    
    header = create_vcf_header(ref_fa, sample_name, filters_from_data)
    
    with pysam.VariantFile(str(output_vcf), 'w', header=header) as vcf_out:
        for row in df.iter_rows(named=True):
            try:
                record = variant_to_vcf_record(row, header, sample_name)
                vcf_out.write(record)
            except Exception as e:
                logging.warning(f"Failed to write variant {row['id']}: {e}")
                continue
    
    logging.info(f"VCF written successfully to {output_vcf}")


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
        help='Output VCF file path (use .vcf.gz for compressed)'
    )
    parser.add_argument(
        '--sample',
        required=True,
        help='Sample name for VCF'
    )
    parser.add_argument(
        '--ref',
        help='Reference FASTA file (for contig headers)'
    )
    
    args = parser.parse_args()
    
    parquet_files = [Path(p) for p in args.parquet]
    output_vcf = Path(args.output)
    ref_fa = Path(args.ref) if args.ref else None
    
    for pf in parquet_files:
        if not pf.exists():
            logging.error(f"Parquet file not found: {pf}")
            sys.exit(1)
    
    if ref_fa and not ref_fa.exists():
        logging.error(f"Reference FASTA not found: {ref_fa}")
        sys.exit(1)
    
    df = load_pav_calls(parquet_files)
    write_vcf(df, output_vcf, args.sample, ref_fa)
    
    if output_vcf.suffix == '.gz':
        logging.info(f"Indexing {output_vcf}")
        pysam.tabix_index(str(output_vcf), preset='vcf', force=True)
        logging.info("Indexing complete")
    
    logging.info("Done!")


if __name__ == '__main__':
    main()
