#!/usr/bin/env python3

import argparse
from pysam import VariantFile


def main():
    parser = argparse.ArgumentParser(description="Revert symbolic alleles to their original representation.")
    parser.add_argument("--input", required=True, help="Path to the annotated VCF (symbolic)")
    parser.add_argument("--output", required=True, help="Path to output VCF (default: stdout)")
    parser.add_argument("--original", required=True, help="Path to the original VCF (non-symbolic)")
    args = parser.parse_args()

    original_data = {}

    with VariantFile(args.original) as original_vcf:
        for record in original_vcf:
            svlen = record.info.get("SVLEN")
            if svlen and isinstance(svlen, (list, tuple)):
                svlen = svlen[0]
            original_data[record.id] = (record.ref, record.alts, svlen)

    vcf_in = VariantFile(args.annotated)
    vcf_out = VariantFile(args.output, 'w', header=vcf_in.header)

    for record in vcf_in:
        if record.id in original_data:
            ref, alts, svlen = original_data[record.id]

            record.ref = ref
            record.alts = alts
            if svlen is not None:
                record.info["SVLEN"] = svlen

            if "BND_ALT" in record.info:
                del record.info["BND_ALT"]

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
