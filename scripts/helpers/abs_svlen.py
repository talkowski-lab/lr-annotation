import sys
from pysam import VariantFile


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile("-", "w", header=vcf_in.header)

    for rec in vcf_in.fetch():
        if "SVLEN" in rec.info:
            if isinstance(rec.info["SVLEN"], tuple):
                rec.info["SVLEN"] = (abs(rec.info["SVLEN"][0]),)
            else:
                rec.info["SVLEN"] = abs(rec.info["SVLEN"])
        vcf_out.write(rec)  # write every record


if __name__ == "__main__":
    main()
