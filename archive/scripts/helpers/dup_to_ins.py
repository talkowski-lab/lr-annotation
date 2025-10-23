import sys
from pysam import VariantFile, FastaFile


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile("-", "w", header=vcf_in.header)

    vcf_in.header.add_line(
        '##INFO=<ID=ORIG_DUP,Number=0,Type=Flag,Description="Record was originally called as DUP">'
    )
    vcf_out.header.add_line(
        '##INFO=<ID=ORIG_DUP,Number=0,Type=Flag,Description="Record was originally called as DUP">'
    )

    for rec in vcf_in.fetch():
        # DUP-->INS
        if rec.info["SVTYPE"] == "DUP":
            rec.info["SVTYPE"] = "INS"
            rec.info.__setitem__("ORIG_DUP", 1)  # set flag
            rec.stop = rec.pos + 1
            if "<DUP>" in rec.alts:
                rec.alts = ("<INS>",)
        vcf_out.write(rec)  # write every variant


if __name__ == "__main__":
    main()
