#!/usr/bin/env python3


import sys
from pysam import VariantFile


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: program.py <annotated_vcf> <original_vcf>\n")
        sys.exit(1)

    annotated_vcf_path = sys.argv[1]
    original_vcf_path = sys.argv[2]

    with VariantFile(annotated_vcf_path) as annotated_vcf, VariantFile(
        original_vcf_path
    ) as original_vcf:
        original_data = {}
        for record in original_vcf.fetch():
            original_data[record.id] = (record.ref, record.alts, record.info.get("SVLEN"))

        vcf_out = VariantFile("-", "w", header=annotated_vcf.header)
        for annotated_record in annotated_vcf.fetch():
            if annotated_record.id in original_data:
                ref, alts, svlen = original_data[annotated_record.id]

                new_record = annotated_record.copy()
                new_record.ref = ref
                new_record.alts = alts
                if svlen is not None:
                    new_record.info["SVLEN"] = svlen

                if "BND_ALT" in new_record.info:
                    del new_record.info["BND_ALT"]

                vcf_out.write(new_record)
            else:
                vcf_out.write(annotated_record)


if __name__ == "__main__":
    main()
