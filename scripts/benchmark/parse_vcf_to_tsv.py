import argparse
import sys
from pysam import VariantFile

def parse_vcf_to_tsv(vcf_path, output_path):
    """
    Parses a VCF file and extracts key information into a TSV file.
    """
    try:
        vcf_in = VariantFile(vcf_path)
    except Exception as e:
        sys.exit(f"Error opening VCF file {vcf_path}: {e}")

    with open(output_path, 'w') as f_out:
        # Write header
        f_out.write("CHROM\tPOS\tREF\tALT\tINFO\n")

        for rec in vcf_in:
            # Reconstruct the INFO field as a string
            info_str = ";".join(f"{key}={','.join(map(str, val)) if isinstance(val, tuple) else str(val)}" for key, val in rec.info.items())
            
            # Write variant information to TSV
            for alt in rec.alts:
                f_out.write(f"{rec.chrom}\t{rec.pos}\t{rec.ref}\t{alt}\t{info_str}\n")

    vcf_in.close()

def main():
    parser = argparse.ArgumentParser(description="Parse a VCF file and output a condensed TSV.")
    parser.add_argument("vcf_path", help="Path to the input VCF file.")
    parser.add_argument("output_path", help="Path to the output TSV file.")
    args = parser.parse_args()

    parse_vcf_to_tsv(args.vcf_path, args.output_path)

if __name__ == "__main__":
    main() 