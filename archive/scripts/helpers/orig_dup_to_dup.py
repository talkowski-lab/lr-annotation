import sys
from pysam import VariantFile, FastaFile


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile("-", "w", header=vcf_in.header)

    for rec in vcf_in.fetch():
        if rec.info["SVTYPE"] == "INS" and "ORIG_DUP" in rec.info:
            rec.info["SVTYPE"] = "DUP"
            if isinstance(rec.info["SVLEN"], tuple):
                svlen = rec.info["SVLEN"][0]
            else:
                svlen = rec.info["SVLEN"]
            rec.stop = rec.pos + svlen
            if "<INS>" in rec.alts:
                rec.alts = ("<DUP>",)
        vcf_out.write(rec)  # write everything


if __name__ == "__main__":
    main()
