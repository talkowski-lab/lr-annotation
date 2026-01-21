#!/usr/bin/env python3


import sys
from pysam import VariantFile


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: program.py vcf svtypes_file\n")
        sys.exit(1)
    
    vcf_path = sys.argv[1]
    svtypes_file = sys.argv[2]
    
    with open(svtypes_file) as f:
        present_svtypes = set(line.strip() for line in f if line.strip())
    
    vcf_in = VariantFile(vcf_path)
    header = vcf_in.header
    
    if len(present_svtypes) > 0:
        header.add_line('##ALT=<ID=N,Description="Baseline reference">')
    
    if "BND" in present_svtypes:
        header.add_line(
            '##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">'
        )
    
    alt_definitions = {
        "DEL": '##ALT=<ID=DEL,Description="Deletion">',
        "INS": '##ALT=<ID=INS,Description="Insertion">',
        "DUP": '##ALT=<ID=DUP,Description="Duplication">',
        "INV": '##ALT=<ID=INV,Description="Inversion">'
    }
    
    for alt_id, alt_line in alt_definitions.items():
        if alt_id in present_svtypes and alt_id not in header.alts:
            header.add_line(alt_line)
    
    vcf_out = VariantFile("-", "w", header=header)

    for rec in vcf_in.fetch():
        if rec.info["SVTYPE"] == "BND":
            rec.info["BND_ALT"] = rec.alts[0]

        rec.alts = ("<%s>" % rec.info["SVTYPE"],)
        rec.ref = "N"

        if "SVLEN" in rec.info:
            if isinstance(rec.info["SVLEN"], tuple):
                rec.info["SVLEN"] = (abs(rec.info["SVLEN"][0]),)
            else:
                rec.info["SVLEN"] = abs(rec.info["SVLEN"])

        if rec.info["SVTYPE"] in ["DEL", "DUP", "INV"]:
            if isinstance(rec.info["SVLEN"], tuple):
                svlen = abs(rec.info["SVLEN"][0])
            else:
                svlen = abs(rec.info["SVLEN"])
            rec.stop = rec.pos + svlen

        vcf_out.write(rec)


if __name__ == "__main__":
    main()
