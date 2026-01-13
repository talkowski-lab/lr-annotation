#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pysam
import sys


def merge_headers(vcf_files):
    base_header = vcf_files[0].header.copy()
    for vcf in vcf_files[1:]:
        for key, info in vcf.header.info.items():
            if key not in base_header.info:
                header_line = f'##INFO=<ID={info.name},Number={info.number},Type={info.type},Description="{info.description}">'
                base_header.add_line(header_line)
    return base_header


def merge_vcf_files(input_files, output_file):
    vcf_files = [pysam.VariantFile(f) for f in input_files]
    merged_header = merge_headers(vcf_files)
    output_vcf = pysam.VariantFile(output_file, "w", header=merged_header)

    iterators = [iter(vcf) for vcf in vcf_files]

    for records in zip(*iterators):
        positions = [(r.chrom, r.pos) for r in records]
        if len(set(positions)) != 1:
            print(f"Error: Records don't match - {positions}")
            sys.exit(1)

        base_record = records[0]
        new_record = merged_header.new_record(
            contig=base_record.chrom,
            start=base_record.start,
            stop=base_record.stop,
            alleles=base_record.alleles,
            id=base_record.id,
            qual=base_record.qual,
            filter=base_record.filter,
        )

        for sample_name in base_record.samples:
            new_record.samples[sample_name].update(base_record.samples[sample_name])

        for record in records:
            for key, value in record.info.items():
                if key not in new_record.info:
                    new_record.info[key] = value

        output_vcf.write(new_record)

    for vcf in vcf_files:
        vcf.close()
    output_vcf.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcfs", nargs="+", help="Input VCF files to merge")
    parser.add_argument("-o", "--output", required=True, help="Output merged VCF file")

    args = parser.parse_args()
    merge_vcf_files(args.input_vcfs, args.output)
    print(f"Merged {len(args.input_vcfs)} VCF files into {args.output}")


if __name__ == "__main__":
    main()
