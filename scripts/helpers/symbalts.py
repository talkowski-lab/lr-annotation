#!/usr/bin/env python3

import argparse
from pysam import VariantFile

NULL_GT = [(0, 0), (None, None), (0, ), (None, ), (None, 0)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
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

    for svtype in present_types:
        if svtype not in header.alts:
            header.add_line(f'##ALT=<ID={svtype},Description="{svtype} variant">')
        if svtype == "BND" and "BND_ALT" not in header.info:
            header.add_line('##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">')

    vcf_out = VariantFile(args.output, 'w', header=header)
    for record in vcf_in:
        if "SVTYPE" in record.info:
            svtype = record.info["SVTYPE"]
            if isinstance(svtype, (list, tuple)):
                svtype = svtype[0]

            svtype = svtype.upper()
            if svtype == "BND":
                record.info["BND_ALT"] = record.alts[0]
            record.ref = "N"
            record.alts = (f"<{svtype}>", )

        if "SVLEN" in record.info:
            svlen = record.info["SVLEN"]
            if isinstance(svlen, (list, tuple)):
                svlen = svlen[0]

            svlen = abs(svlen)
            record.stop = record.pos + svlen

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
