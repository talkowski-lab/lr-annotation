#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(
        description="Find exact matches between two VCF files."
    )
    parser.add_argument("vcf_eval", help="Evaluation VCF file.")
    parser.add_argument("vcf_truth", help="Truth VCF file.")
    parser.add_argument("prefix", help="Prefix for output files.")
    args = parser.parse_args()

    truth_variants = {}
    with pysam.VariantFile(args.vcf_truth) as vcf_in:
        for record in vcf_in:
            key = (record.chrom, record.pos, record.ref, record.alts)
            truth_variants[key] = record.id

    vcf_in = pysam.VariantFile(args.vcf_eval)

    vcf_in.header.info.add(
        "gnomAD_V4_match", "1", "String", "Matching status against gnomAD v4."
    )
    vcf_in.header.info.add(
        "gnomAD_V4_match_ID", "1", "String", "Matching variant ID from gnomAD v4."
    )

    matched_tmp_path = f"{args.prefix}.exact_matched.tmp.vcf"
    unmatched_tmp_path = f"{args.prefix}.unmatched.tmp.vcf"
    matched_out_path = f"{args.prefix}.exact_matched.vcf.gz"
    unmatched_out_path = f"{args.prefix}.unmatched.vcf.gz"

    with open(matched_tmp_path, "w") as matched_out, open(
        unmatched_tmp_path, "w"
    ) as unmatched_out:
        matched_out.write(str(vcf_in.header))
        unmatched_out.write(str(vcf_in.header))

        for record in vcf_in:
            key = (record.chrom, record.pos, record.ref, record.alts)
            if key in truth_variants:
                record.info["gnomAD_V4_match"] = "EXACT"
                record.info["gnomAD_V4_match_ID"] = truth_variants[key]
                matched_out.write(str(record))
            else:
                unmatched_out.write(str(record))

    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", matched_out_path, matched_tmp_path],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", "-f", matched_out_path], check=True)
    os.remove(matched_tmp_path)

    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", unmatched_out_path, unmatched_tmp_path],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", "-f", unmatched_out_path], check=True)
    os.remove(unmatched_tmp_path)


if __name__ == "__main__":
    main()
