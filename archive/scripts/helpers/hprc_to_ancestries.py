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
        "--ancestries_reference",
        required=False,
        help="Path to the ancestries reference TSV file (fallback mapping)"
    )
    parser.add_argument(
        "--exclude_samples",
        default="",
        help="Comma-separated list of sample IDs to exclude"
    )
    parser.add_argument(
        "--swap_pop",
        action="append",
        help="Swap population names in output (format: old,new). Can be specified multiple times."
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output file (default: stdout)"
    )
    
    args = parser.parse_args()
    
    swap_map = {}
    if args.swap_pop:
        for swap in args.swap_pop:
            parts = swap.split(",")
            if len(parts) == 2:
                old_pop, new_pop = parts
                swap_map[old_pop.lower()] = new_pop.lower()
            else:
                print(f"Warning: Invalid swap format '{swap}', expected 'old,new'", file=sys.stderr)
    
    exclude_set = set(args.exclude_samples.split(",")) if args.exclude_samples else set()
    exclude_set.discard("")
    
    ancestry_map = {}
    with open(args.ancestries, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            ethnicity = row["ethnicity"].lower()
            ancestry_map[sample] = ethnicity
    
    pop_to_ancestry = {}
    if args.ancestries_reference:
        with open(args.ancestries_reference, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) >= 2:
                    ancestry = row[0].lower()
                    population_abbrev = row[1]
                    pop_to_ancestry[population_abbrev] = ancestry
    
    output_file = sys.stdout if args.output == "-" else open(args.output, "w")
    
    with open(args.metadata, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            sample_id = row["sample_id"]
            
            if sample_id in exclude_set:
                continue
            
            if sample_id in ancestry_map:
                ethnicity = ancestry_map[sample_id]
                ethnicity = swap_map.get(ethnicity.lower(), ethnicity)
                output_file.write(f"{sample_id}\t{ethnicity}\n")
            else:
                population_abbrev = row["population_abbreviation"]
                if population_abbrev in pop_to_ancestry:
                    ethnicity = pop_to_ancestry[population_abbrev]
                    ethnicity = swap_map.get(ethnicity.lower(), ethnicity)
                    output_file.write(f"{sample_id}\t{ethnicity}\n")
                else:
                    print(f"{sample_id}")
    
    if output_file != sys.stdout:
        output_file.close()


if __name__ == "__main__":
    main()
