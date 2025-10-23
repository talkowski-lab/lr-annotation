import sys
from pysam import VariantFile, FastaFile


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: program.py vcf")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile("-", "w", header=vcf_in.header)

    cur_pos = -1
    cur_prev = []
    for rec in vcf_in.fetch():
        # INS-->DUP
        if rec.info["SVTYPE"] not in ["INS", "DUP"]:
            vcf_out.write(rec)
            continue
        if not cur_prev:
            cur_pos = rec.pos
            cur_prev = [rec]
            continue
        # match on chrom, pos, and svlen and mismatch on svtype (can't match on ref, pos after standardizing/annotation!)
        if cur_pos == rec.pos:
            for prev in cur_prev:
                if (prev.info["SVLEN"] == rec.info["SVLEN"]) and prev.info[
                    "SVTYPE"
                ] != rec.info["SVTYPE"]:
                    if prev.info["SVTYPE"] == "INS":
                        keep = prev
                        discard = rec
                    else:
                        keep = rec
                        discard = prev

                    # now collapse the annotations
                    for anno in discard.info.keys():
                        if not anno.startswith("PREDICTED"):
                            continue  # only collapse the predicted consequence annotations
                        if anno == "PREDICTED_INTERGENIC":
                            continue  # this is a boolean, not interesting to transfer
                        if anno not in keep.info:
                            keep.info[anno] = discard.info[anno]
                        else:
                            keep.info[anno] = tuple(
                                set(keep.info[anno] + discard.info[anno])
                            )  # merge and keep unique entries
                    vcf_out.write(keep)
            cur_prev.append(rec)
        else:  # new pos
            cur_pos = rec.pos
            cur_prev = [rec]


if __name__ == "__main__":
    main()
