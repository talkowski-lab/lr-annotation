#!/usr/bin/env python3


import sys
from Bio import SeqIO
from pysam import VariantFile

# intersection output:
# 0:chrom 1:range_start 2:range_end 3:position(s) 4:REF(s) 5:SVLEN(s) (<<
# from target VCF <<) 6:chrom 7:pos-1 8:pos 9:mei_len 10:sample(s) (<<
# from MEI VCF <<)


def main():
    if len(sys.argv) != 5:
        sys.stderr.write(
            "Usage: program.py intersect_output target_VCF INS_fa ME_type\n"
        )
        sys.exit(1)

    intersect_output = open(sys.argv[1], "r")
    target_VCF = VariantFile(sys.argv[2])
    INS_fa = SeqIO.index(sys.argv[3], "fasta")
    ME_type = sys.argv[4]

    if ME_type not in ["SVA", "ALU", "Alu", "L1", "LINE1", "LINE", "HERVK"]:
        sys.stderr.write("ME type %s not recognized\n" % ME_type)
        sys.exit(1)

    samples = list(target_VCF.header.samples)

    for line in intersect_output:
        line = line.strip().split("\t")
        chrom = line[0]
        MEI_length = int(line[9])
        MEI_samples = line[10].rstrip(",").split(",")

        SVs = []
        SV_positions = line[3].split(",")
        SV_refs = line[4].split(",")
        SV_lengths = line[5].split(",")

        for ref, pos, length in zip(SV_refs, SV_positions, SV_lengths):
            SVs.append(
                {
                    "ref": ref.split("_")[0],
                    "pos": int(pos),
                    "length": int(length),
                    "id": "%s:%s;%s" % (chrom, pos, ref),
                }
            )

        for SV in SVs:
            # match lengths on 80% reciprocal size percentage but ignore if SVA
            ratio = min(SV["length"], MEI_length) / max(SV["length"], MEI_length)
            if ME_type == "SVA" or ratio >= 0.8:
                for rec in target_VCF.fetch(
                    chrom, SV["pos"] - 1, SV["pos"]
                ):  # tested that this is the right way to fetch one position
                    # skip if the record already has a ME_TYPE annotation
                    if "ME_TYPE" in rec.info:
                        continue

                    # match on svlen - 1 less as alt length includes the ref nuc
                    if isinstance(rec.info["SVLEN"], tuple):
                        svlen = abs(rec.info["SVLEN"][0])
                    else:
                        svlen = abs(rec.info["SVLEN"])
                    if svlen != SV["length"] - 1:
                        continue

                    # ensure that the SVTYPE is an insertion
                    if isinstance(rec.info["SVTYPE"], tuple):
                        svtype = rec.info["SVTYPE"][0]
                    else:
                        svtype = rec.info["SVTYPE"]
                    if svtype != "INS":
                        continue

                    # match on ref
                    if rec.ref != SV["ref"]:
                        continue

                    # match on alt
                    if rec.alts[0] != str(INS_fa[SV["id"]].seq):
                        continue

                    # check that there's enough sample overlap
                    GTs = [rec.samples[sample]["GT"] for sample in samples]
                    SV_samples = []
                    for GT, sample in zip(GTs, samples):
                        if 1 in GT:
                            SV_samples.append(sample)
                    shared_samples = list(set(SV_samples) & set(MEI_samples))
                    if (
                        shared_samples
                    ):  # at least one shared sample between MEI and SV calls
                        # write annotation file:
                        sys.stdout.write(
                            "%s\t%s\t%s\t%s\t%s\t%d\n"
                            % (
                                chrom,
                                SV["pos"],
                                rec.ref,
                                rec.alts[0],
                                ME_type,
                                MEI_length,
                            )
                        )
                    else:
                        sys.stderr.write("%s\n%s\n#####\n" % (MEI_samples, SV_samples))


if __name__ == "__main__":
    main()
