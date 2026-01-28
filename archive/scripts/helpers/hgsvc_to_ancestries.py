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
    parser = argparse.ArgumentParser(description="Create an ancestry file from mapping and VCF samples")
    parser.add_argument("--vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--ancestries", required=True, help="Path to ancestries TSV")
    parser.add_argument("--swap_samples", help="Comma-separated pair: old_id,new_id")
    parser.add_argument("--swap_pop", action="append", help="Swap population names in output (format: old,new)")
    parser.add_argument("--output", default="-", help="Output file (default: stdout)")
    
    args = parser.parse_args()

    pop_swap_map = {}
    if args.swap_pop:
        for swap in args.swap_pop:
            parts = swap.split(",")
            if len(parts) == 2:
                old_pop, new_pop = parts
                pop_swap_map[old_pop.lower()] = new_pop.lower()
            else:
                print(f"Warning: Invalid swap format '{swap}', expected 'old,new'", file=sys.stderr)
        
    swap_old, swap_new = None, None
    if args.swap_samples:
        parts = args.swap_samples.split(',')
        if len(parts) == 2:
            swap_old, swap_new = parts
            
    vcf_samples = get_vcf_samples(args.vcf)
    
    ancestry_map = {}
    with open(args.ancestries, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get("sample")
            ethnicity = row.get("ethnicity")
            if sample and ethnicity:
                ancestry_map[sample] = ethnicity.lower()
    
    output_file = sys.stdout if args.output == "-" else open(args.output, "w")
    
    for sample_id in vcf_samples:
        if sample_id in ancestry_map:
            active_id = swap_new if sample_id == swap_old else sample_id
            ethnicity = ancestry_map[sample_id]
            ethnicity = pop_swap_map.get(ethnicity.lower(), ethnicity)
            output_file.write(f"{active_id}\t{ethnicity}\n")
        else:
            print(f"Sample {sample_id} not found in ancestries - skipping.")
    
    if output_file != sys.stdout:
        output_file.close()

if __name__ == "__main__":
    main()