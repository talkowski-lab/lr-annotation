#!/usr/bin/env python3

import argparse
import pysam
import bgzip

def main():
    parser = argparse.ArgumentParser(description="Find exact matches between two VCF files.")
    parser.add_argument("vcf_eval", help="Evaluation VCF file.")
    parser.add_argument("vcf_truth", help="Truth VCF file.")
    parser.add_argument("prefix", help="Prefix for output files.")
    args = parser.parse_args()

    truth_variants = set()
    with pysam.VariantFile(args.vcf_truth) as vcf_in:
        for record in vcf_in:
            key = (record.chrom, record.pos, record.ref, record.alts)
            truth_variants.add(key)

    vcf_in = pysam.VariantFile(args.vcf_eval)
    
    # Add header for the new INFO field
    vcf_in.header.info.add('gnomAD_V4_match', '1', 'String', 'Matching status against gnomAD v4.')

    matched_out_path = f"{args.prefix}.exact_matched.vcf.gz"
    unmatched_out_path = f"{args.prefix}.unmatched.vcf.gz"

    with bgzip.BGZFile(matched_out_path, "w") as matched_out, bgzip.BGZFile(unmatched_out_path, "w") as unmatched_out:
        matched_out.write(str(vcf_in.header).encode())
        unmatched_out.write(str(vcf_in.header).encode())
        
        for record in vcf_in:
            key = (record.chrom, record.pos, record.ref, record.alts)
            if key in truth_variants:
                record.info['gnomAD_V4_match'] = 'EXACT'
                matched_out.write(str(record).encode())
            else:
                unmatched_out.write(str(record).encode())

    pysam.tabix_index(matched_out_path, preset="vcf")
    pysam.tabix_index(unmatched_out_path, preset="vcf")

if __name__ == "__main__":
    main() 