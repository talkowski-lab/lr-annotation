#!/usr/bin/env python3


import sys
from pysam import VariantFile


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile("-", "w", header=vcf_in.header)

    vcf_in.header.add_line(
        '##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">'
    )
    vcf_out.header.add_line(
        '##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">'
    )
    
    alt_definitions = {
        "DEL": '##ALT=<ID=DEL,Description="Deletion">',
        "DUP": '##ALT=<ID=DUP,Description="Duplication">',
        "INV": '##ALT=<ID=INV,Description="Inversion">',
        "INS": '##ALT=<ID=INS,Description="Insertion">'
    }
    
    for alt_id, alt_line in alt_definitions.items():
        if alt_id not in vcf_out.header.alts:
            vcf_out.header.add_line(alt_line)

    for rec in vcf_in.fetch():
        if rec.info["SVTYPE"] == "BND":
            rec.info["BND_ALT"] = rec.alts[0]

        rec.alts = ("<%s>" % rec.info["SVTYPE"],)
        rec.ref = "N"

        if rec.info["SVTYPE"] in ["DEL", "DUP", "INV"]:
            if isinstance(rec.info["SVLEN"], tuple):
                svlen = abs(rec.info["SVLEN"][0])
            else:
                svlen = abs(rec.info["SVLEN"])
            rec.stop = rec.pos + svlen

        vcf_out.write(rec)


if __name__ == "__main__":
    main()
