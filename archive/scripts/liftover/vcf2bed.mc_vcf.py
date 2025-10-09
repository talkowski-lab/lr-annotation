#!/usr/bin/env python3
import sys
import gzip

def vcf_to_bed(vcf_file, bed_file):
    # open input (support gzipped or plain vcf)
    if vcf_file.endswith(".gz"):
        infile = gzip.open(vcf_file, "rt")
    else:
        infile = open(vcf_file, "r")
    
    with open(bed_file, "w") as out:
        for line in infile:
            if line.startswith("#"):
                continue  # skip header lines
            
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1]) - 1  # BED is 0-based
            ref = fields[3]
            alt = fields[4]
            end = pos + len(ref)      # END = POS + len(REF)

            # write BED: chrom, start, end, ref, alt
            out.write(f"{chrom}\t{pos}\t{end}\t{ref}\t{alt}\n")

    infile.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vcf_to_bed.py input.vcf[.gz] output.bed")
        sys.exit(1)

    vcf_file = sys.argv[1]
    bed_file = sys.argv[2]
    vcf_to_bed(vcf_file, bed_file)
