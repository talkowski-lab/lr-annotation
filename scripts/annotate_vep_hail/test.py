from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os
import json
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='vcf_file', help='Input VCF file')
parser.add_argument('-o', dest='vep_annotated_vcf_name', help='Output filename')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--project-id', dest='project_id', help='Google Project ID')

args = parser.parse_args()
vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores
mem = int(np.floor(float(args.mem)))
build = args.build
gcp_project = args.project_id

hl.init(
    min_block_size=128, 
    spark_conf = {
        "spark.executor.cores": cores, 
        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
        "spark.driver.cores": cores,
        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
    }, 
    tmp_dir="tmp", 
    local_tmpdir="tmp",
)

print(f"Reading VCF: {vcf_file}")
header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

print("Running VEP...")
mt = hl.vep(mt, config='vep_config.json', csq=False, tolerate_parse_error=True)

print("Selecting the most severe consequence from the VEP struct...")
most_severe_csq = hl.vep_transcript_consequences(mt.vep).most_severe()

print("Flattening VEP struct to match gnomAD format...")
vep_struct = hl.struct(
    transcript_consequences=hl.array([
        hl.struct(
            feature_type=most_severe_csq.biotype,
            feature_id=most_severe_csq.transcript_id,
            consequence_terms=most_severe_csq.consequence_terms,
            impact=most_severe_csq.impact,
            gene_id=most_severe_csq.gene_id,
            gene_symbol=most_severe_csq.gene_symbol,
            canonical=most_severe_csq.canonical,
            polyphen_prediction=most_severe_csq.polyphen_prediction,
            sift_prediction=most_severe_csq.sift_prediction,
            hgvsc=most_severe_csq.hgvsc,
            hgvsp=most_severe_csq.hgvsp
        )
    ]),
    most_severe_consequence=mt.vep.most_severe_consequence,
    variant_class=mt.vep.variant_class,
)

# STEP 4: FINAL STEP BEFORE EXPORT - Serialize the curated struct into a JSON string in the INFO field.
print("Annotating final JSON string to INFO field...")
mt = mt.annotate_rows(info=mt.info.annotate(vep=hl.json(vep_struct)))

# STEP 5: Clean up by dropping the large, now-redundant temporary fields.
mt = mt.drop('vep', 'vep_proc_id')

# STEP 6: Dynamically create the VCF header description string.
print("Updating VCF header...")
vep_field_description = " | ".join(
    [f"{f}: {t}" for f, t in vep_struct.transcript_consequences.dtype.element_type.items()]
)
vep_description = (
    "Consequence annotations from Ensembl VEP. "
    "Only the most severe consequence across all transcripts is included. "
    f"Format: {vep_field_description}"
)
header['info']['vep'] = {'Description': vep_description, 'Number': '.', 'Type': 'String'}

# STEP 7: Export the final VCF.
print(f"Exporting annotated VCF to: {vep_annotated_vcf_name}")
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)

print("VEP annotation complete.")