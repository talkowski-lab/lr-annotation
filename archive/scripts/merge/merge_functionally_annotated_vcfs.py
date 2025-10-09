#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Merge VEP and SVAnnotate VCF files based on variant size.
"""

import argparse
import pysam
import sys


def get_variant_length(record):
    svlen = record.info['SVLEN']
    if isinstance(svlen, (list, tuple)):
        svlen = svlen[0]
    return abs(svlen)


def merge_headers(vep_vcf, svannotate_vcf):
    base_header = vep_vcf.header.copy()
    
    vep_info_fields = set(vep_vcf.header.info.keys())
    svannotate_info_fields = set(svannotate_vcf.header.info.keys())
    
    # Add SVAnnotate exclusive fields
    for key, info in svannotate_vcf.header.info.items():
        if key not in base_header.info:
            description = info.description
            if description.endswith('.'):
                description = description[:-1]
            if not description.endswith('(only applies to variants ≥50 bp)'):
                description = f'{description} (only applies to variants ≥50 bp)'
            header_line = f'##INFO=<ID={info.name},Number={info.number},Type={info.type},Description="{description}">'
            base_header.add_line(header_line)
    
    # Modify VEP exclusive fields by creating new header lines
    for key, info in vep_vcf.header.info.items():
        if key not in svannotate_info_fields:
            description = info.description
            if description.endswith('.'):
                description = description[:-1]
            if not description.endswith('(only applies to variants <50bp)'):
                description = f'{description} (only applies to variants <50bp)'
            header_line = f'##INFO=<ID={info.name},Number={info.number},Type={info.type},Description="{description}">'
            base_header.add_line(header_line)
    
    return base_header


def merge_vcf_files(vep_file, svannotate_file, output_file):
    vep_vcf = pysam.VariantFile(vep_file)
    svannotate_vcf = pysam.VariantFile(svannotate_file)
    
    merged_header = merge_headers(vep_vcf, svannotate_vcf)
    output_vcf = pysam.VariantFile(output_file, 'w', header=merged_header)
    
    record_count = 0
    small_variant_count = 0
    large_variant_count = 0
    
    for vep_record, svannotate_record in zip(vep_vcf, svannotate_vcf):
        record_count += 1
        vep_pos = (vep_record.chrom, vep_record.pos)
        svannotate_pos = (svannotate_record.chrom, svannotate_record.pos)
        if vep_pos != svannotate_pos:
            print(f"Error: Records don't match at position {record_count}")
            print(f"VEP: {vep_pos}")
            print(f"SVAnnotate: {svannotate_pos}")
            sys.exit(1)
        
        variant_length = get_variant_length(svannotate_record)
        if variant_length >= 50:
            source_record = svannotate_record
            large_variant_count += 1
        else:
            source_record = vep_record
            small_variant_count += 1
        
        new_record = merged_header.new_record(
            contig=source_record.chrom,
            start=source_record.start,
            stop=source_record.stop,
            alleles=source_record.alleles,
            id=source_record.id,
            qual=source_record.qual,
            filter=source_record.filter
        )
        
        for key, value in source_record.info.items():
            new_record.info[key] = value
        
        if hasattr(source_record, 'samples') and source_record.samples:
            for sample_name in source_record.samples:
                new_record.samples[sample_name].update(source_record.samples[sample_name])
        
        output_vcf.write(new_record)
    
    vep_vcf.close()
    svannotate_vcf.close()
    output_vcf.close()
    
    print(f"Merged {record_count} records:")
    print(f"  Small variants (<50bp): {small_variant_count} (from VEP)")
    print(f"  Large variants (≥50bp): {large_variant_count} (from SVAnnotate)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--vep-vcf', required=True,
                       help='VEP-annotated VCF file (for variants <50bp)')
    parser.add_argument('--svannotate-vcf', required=True,
                       help='SVAnnotate-annotated VCF file (for variants ≥50bp)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output merged VCF file')
    
    args = parser.parse_args()
    
    merge_vcf_files(args.vep_vcf, args.svannotate_vcf, args.output)


if __name__ == '__main__':
    main()
