#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <input_csv_file> <output_ped_file>"
    exit 1
fi

INPUT_CSV="$1"
OUTPUT_PED="$2"

export PYTHONIOENCODING=latin-1

python3 -c '
import sys, csv

reader = csv.reader(sys.stdin, delimiter="\t")

next(reader)

for row in reader:
    fam_id = row[0]
    ind_id = row[0]
    pat_id = "0"
    mat_id = "0"
    sex = "1" if row[1] == "male" else "2" if row[1] == "female" else "0"
    pheno = "0"
    print("\t".join([fam_id, ind_id, pat_id, mat_id, sex, pheno]))

' < "$INPUT_CSV" > "$OUTPUT_PED"

echo "PED file saved to $OUTPUT_PED"
