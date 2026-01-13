#!/usr/bin/env python3
import sys
import gzip


def load_bed(bed_file):
    bed_dict = {}
    with gzip.open(bed_file, "rt") if bed_file.endswith(".gz") else open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, ref, alt = line.strip().split()[:5]
            key = (chrom, ref, alt)
            bed_dict[key] = (int(start), int(end))
    return bed_dict


def update_vcf(vcf_file, bed_dict, output_file):
    with gzip.open(vcf_file, "rt") as f_in, gzip.open(output_file, "wt") as f_out:
        for line in f_in:
            if line.startswith("#"):
                f_out.write(line)
                continue

            fields = line.strip().split("\t")
            chrom, pos, vid, ref, alt = fields[:5]
            key = (chrom, ref, alt)

            if key in bed_dict:
                start, end = bed_dict[key]
                fields[1] = str(start + 1)  # BED is 0-based, VCF is 1-based

            f_out.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python bed_vcf_update.py input.bed(.gz) input.vcf.gz output.vcf.gz"
        )
        sys.exit(1)

    bed_file, vcf_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    bed_dict = load_bed(bed_file)
    update_vcf(vcf_file, bed_dict, output_file)
