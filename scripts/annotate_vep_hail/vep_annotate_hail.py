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
mt = mt.annotate_rows(info = mt.info.annotate(vep=hl.json(mt.vep)))
mt = mt.drop('vep', 'vep_proc_id')

print("Selecting most severe consequence...")
most_severe_csq = hl.vep_transcript_consequences(mt.vep).most_severe()

print("Flattening structure...")
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

print("Annotating JSON string onto main matrix table...")
mt = mt.annotate_rows(info=mt.info.annotate(vep=hl.json(vep_struct)))
mt = mt.drop('vep', 'vep_proc_id')

print("Updating VCF header...")
vep_format_string = "|".join(
    [f"{f}" for f, t in vep_struct.transcript_consequences.dtype.element_type.items()]
)
header['info']['vep'] = {'Description': f'Consequence annotations from Ensembl VEP. Format: {vep_format_string}', 'Number': '.', 'Type': 'String'}

print("Exporting annotated VCF...")
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)

print("VEP annotation complete.")