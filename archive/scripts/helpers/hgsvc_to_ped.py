#!/usr/bin/env python3

import argparse
import csv
import sys
import gzip

def get_vcf_samples(vcf_path):
    f = gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')
    for line in f:
        if line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
            f.close()
            return samples
    f.close()
    return []

def main():
    parser = argparse.ArgumentParser(description="Create a PED file from HGSVC metadata and VCF samples")
    parser.add_argument("--vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--metadata", required=True, help="Path to HGSVC metadata TSV")
    parser.add_argument("--swap_samples", help="Comma-separated pair: old_id,new_id")
    parser.add_argument("--output", default="-", help="Output file (default: stdout)")
    
    args = parser.parse_args()
    
    swap_old, swap_new = None, None
    if args.swap_samples:
        parts = args.swap_samples.split(',')
        if len(parts) == 2:
            swap_old, swap_new = parts
    
    vcf_samples = get_vcf_samples(args.vcf)
    
    metadata_map = {}
    with open(args.metadata, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("Sample name")
            if name:
                metadata_map[name] = row

    output_file = sys.stdout if args.output == "-" else open(args.output, "w")
    
    for sample_id in vcf_samples:
        if sample_id not in metadata_map:
            print(f"Sample {sample_id} not found in metadata - skipping.")
            continue
        
        row = metadata_map[sample_id]
        
        active_id = swap_new if sample_id == swap_old else sample_id
        
        sex_str = row.get("Sex", "").lower()
        if sex_str == "male":
            sex_code = "1"
        elif sex_str == "female":
            sex_code = "2"
        else:
            sex_code = "0"
        
        output_file.write(f"{active_id}\t{active_id}\t0\t0\t{sex_code}\t0\n")
            
    if output_file != sys.stdout:
        output_file.close()

if __name__ == "__main__":
    main()