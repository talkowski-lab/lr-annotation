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

reader = csv.reader(sys.stdin)

next(reader)

for row in reader:
    fam_id = row[5] if row[5] else row[0]
    ind_id = row[0]
    pat_id = row[6] if row[6] else "0"
    mat_id = row[7] if row[7] else "0"
    sex_str = row[8]
    
    sex = "1" if sex_str == "male" else "2" if sex_str == "female" else "0"
    
    pheno = "0"
    
    print("\t".join([fam_id, ind_id, pat_id, mat_id, sex, pheno]))

' < "$INPUT_CSV" > "$OUTPUT_PED"

echo "PED file saved to $OUTPUT_PED"
