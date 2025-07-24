import argparse
import sys
from pysam import VariantFile

def get_vep_format_from_header(vcf_path):
    """
    Parses a VCF header to get the VEP/CSQ format string.
    """
    try:
        vcf_in = VariantFile(vcf_path)
    except Exception as e:
        sys.exit(f"Error opening VCF file {vcf_path}: {e}")

    vep_key = None
    vep_format = None

    for header_rec in vcf_in.header.records:
        if header_rec.key == 'INFO' and ('ID' in header_rec and (header_rec['ID'] == 'vep' or header_rec['ID'] == 'CSQ')):
            vep_key = header_rec['ID']
            description = header_rec['Description']
            vep_format = description.split('Format: ')[-1].strip('">')
            break

    vcf_in.close()
    
    if vep_key and vep_format:
        print(f"{vep_key}|{vep_format}")
    else:
        sys.exit("VEP/CSQ format not found in VCF header.")

def main():
    parser = argparse.ArgumentParser(description="Get VEP/CSQ format from VCF header.")
    parser.add_argument("vcf_path", help="Path to the input VCF file.")
    args = parser.parse_args()
    get_vep_format_from_header(args.vcf_path)

if __name__ == "__main__":
    main() 