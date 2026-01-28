#!/usr/bin/env python3

import argparse
from pysam import VariantFile

NULL_GT = [(0, 0), (None, None), (0, ), (None, ), (None, 0)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--svtypes", required=True)
    args = parser.parse_args()

    with open(args.svtypes) as f:
        present_svtypes = set(line.strip() for line in f if line.strip())

    vcf_in = VariantFile(args.input)
    header = vcf_in.header

    if 'END' not in vcf_in.header.info:
        header.add_line('##INFO=<ID=END,Number=.,Type=Integer,Description="End position of the variant">')

    if present_svtypes:
        header.add_line('##ALT=<ID=N,Description="Baseline reference">')

    if "BND" in present_svtypes:
        header.add_line('##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">')

    for svtype in present_svtypes:
        alt_id = svtype
        if alt_id not in header.alts:
            header.add_line(f'##ALT=<ID={alt_id},Description="{alt_id} variant">')

    vcf_out = VariantFile(args.output, 'w', header=header)

    for record in vcf_in:
        for sample in record.samples.values():
            if sample['GT'] in NULL_GT:
                continue

            new_gt = []
            for allele in sample['GT']:
                if allele is None:
                    new_gt.append(None)
                elif allele > 0:
                    new_gt.append(1)
                else:
                    new_gt.append(0)
            sample['GT'] = tuple(new_gt)

        if "SVTYPE" in record.info:
            svtype = record.info["SVTYPE"]
            if isinstance(svtype, (list, tuple)):
                svtype = svtype[0]

            if svtype == "BND":
                record.info["BND_ALT"] = record.alts[0]

            record.ref = "N"
            record.alts = (f"<{svtype}>", )

        if "SVLEN" in record.info:
            svlen = record.info["SVLEN"]
            if isinstance(svlen, (list, tuple)):
                svlen = abs(svlen[0])
            else:
                svlen = abs(svlen)

            record.info["SVLEN"] = svlen
            record.stop = record.pos + svlen

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
