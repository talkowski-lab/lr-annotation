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
        # "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
        # "spark.hadoop.fs.gs.requester.pays.buckets": "hail-datasets-us-central1",
        # "spark.hadoop.fs.gs.requester.pays.project.id": gcp_project,
    }, 
    tmp_dir="tmp", 
    local_tmpdir="tmp",
)

header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

mt = hl.vep(mt, config='vep_config.json', csq=False, tolerate_parse_error=True)
mt = mt.annotate_rows(info = mt.info.annotate(vep=hl.json(mt.vep)))
mt = mt.drop('vep', 'vep_proc_id')

with open('vep_config.json', 'r') as f:
    vep_config = json.load(f)
header['info']['vep'] = {'Description': f"Consequence annotations from Ensembl VEP. Format: {vep_config['vep_json_schema']}", 'Number': '.', 'Type': 'String'}

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)