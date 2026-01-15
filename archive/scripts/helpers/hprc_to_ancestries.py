#!/usr/bin/env python3

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Create an ancestry file from HPRC metadata and ancestry mapping"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to the metadata CSV file"
    )
    parser.add_argument(
        "--ancestries",
        required=True,
        help="Path to the ancestries TSV file"
    )
    parser.add_argument(
        "--exclude_samples",
        default="",
        help="Comma-separated list of sample IDs to exclude"
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output file (default: stdout)"
    )
    
    args = parser.parse_args()
    
    exclude_set = set(args.exclude_samples.split(",")) if args.exclude_samples else set()
    exclude_set.discard("")
    
    ancestry_map = {}
    with open(args.ancestries, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            ethnicity = row["ethnicity"].lower()
            ancestry_map[sample] = ethnicity
    
    output_file = sys.stdout if args.output == "-" else open(args.output, "w")
    
    try:
        with open(args.metadata, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                sample_id = row["sample_id"]
                
                if sample_id in exclude_set:
                    continue
                
                if sample_id in ancestry_map:
                    ethnicity = ancestry_map[sample_id]
                    output_file.write(f"{sample_id}\t{ethnicity}\n")
                else:
                    print(f"{sample_id}")
    
    finally:
        if output_file != sys.stdout:
            output_file.close()


if __name__ == "__main__":
    main()
