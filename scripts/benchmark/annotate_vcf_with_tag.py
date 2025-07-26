import argparse
import pysam
import subprocess
import os

def main():
    parser = argparse.ArgumentParser(description="Add a tag to a VCF file.")
    parser.add_argument("vcf_in", help="Input VCF file.")
    parser.add_argument("vcf_out", help="Output VCF file.")
    parser.add_argument("tag_name", help="Name of the INFO tag to add.")
    parser.add_argument("tag_value", help="Value to assign to the tag.")
    args = parser.parse_args()

    vcf_in = pysam.VariantFile(args.vcf_in)
    vcf_in.header.info.add(args.tag_name, '1', 'String', 'Matching status against gnomAD v4.')

    tmp_out = f"{args.vcf_out}.tmp"

    with open(tmp_out, "w") as vcf_out_handle:
        vcf_out_handle.write(str(vcf_in.header))
        for record in vcf_in:
            record.info[args.tag_name] = args.tag_value
            vcf_out_handle.write(str(record))
    
    subprocess.run(["bcftools", "view", "-Oz", "-o", args.vcf_out, tmp_out], check=True)
    subprocess.run(["tabix", "-p", "vcf", args.vcf_out], check=True)
    os.remove(tmp_out)

if __name__ == "__main__":
    main() 