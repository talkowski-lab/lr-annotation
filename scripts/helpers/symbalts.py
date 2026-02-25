#!/usr/bin/env python3

import argparse
from pysam import VariantFile

NULL_GT = [(0, 0), (None, None), (0, ), (None, ), (None, 0)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--type_field", required=True)
    parser.add_argument("--length_field", required=True)
    parser.add_argument("--types", required=True)
    args = parser.parse_args()

    with open(args.types) as f:
        present_types = set(line.strip().upper() for line in f if line.strip())

    vcf_in = VariantFile(args.input)
    header = vcf_in.header

    if 'END' not in vcf_in.header.info:
        header.add_line('##INFO=<ID=END,Number=.,Type=Integer,Description="End position of the variant">')

    if present_types:
        header.add_line('##ALT=<ID=N,Description="Baseline reference">')

    for allele_type in present_types:
        if allele_type not in header.alts:
            header.add_line(f'##ALT=<ID={allele_type},Description="{allele_type} variant">')
        if allele_type == "BND" and "BND_ALT" not in header.info:
            header.add_line('##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">')

    vcf_out = VariantFile(args.output, 'w', header=header)
    for record in vcf_in:
        # Parse allele type to set REF/ALT
        allele_type = record.info[args.type_field]
        if isinstance(allele_type, (list, tuple)):
            allele_type = allele_type[0]

        allele_type = allele_type.upper()
        if allele_type == "BND":
            record.info["BND_ALT"] = record.alts[0]
        record.ref = "N"
        record.alts = (f"<{allele_type}>", )

        # Parse allele length to set END
        allele_length = record.info[args.length_field]
        if isinstance(allele_length, (list, tuple)):
            allele_length = allele_length[0]

        allele_length = abs(allele_length)
        record.stop = record.pos + allele_length

        # Output record
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
