#!/usr/bin/env python3

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Create a PED file from HPRC metadata CSV"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to the metadata CSV file"
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
    
    output_file = sys.stdout if args.output == "-" else open(args.output, "w")
    
    try:
        with open(args.metadata, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                sample_id = row["sample_id"]
                
                if sample_id in exclude_set:
                    continue
                
                family_id = row["family_id"]
                paternal_id = row["paternal_id"]
                maternal_id = row["maternal_id"]
                sex = row["sex"]
                
                if sex.lower() == "male":
                    sex_code = "1"
                elif sex.lower() == "female":
                    sex_code = "2"
                else:
                    sex_code = "0"
                
                phenotype = "0"
                
                output_file.write(f"{family_id}\t{sample_id}\t{paternal_id}\t{maternal_id}\t{sex_code}\t{phenotype}\n")
    
    finally:
        if output_file != sys.stdout:
            output_file.close()


if __name__ == "__main__":
    main()
