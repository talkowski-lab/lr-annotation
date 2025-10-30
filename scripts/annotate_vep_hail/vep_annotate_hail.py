#!/usr/bin/env python3


import hail as hl
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="Parse arguments")
parser.add_argument("-i", dest="vcf_file", help="Input VCF file")
parser.add_argument("-o", dest="vep_annotated_vcf_name", help="Output filename")
parser.add_argument("--cores", dest="cores", help="CPU cores")
parser.add_argument("--mem", dest="mem", help="Memory")
parser.add_argument("--build", dest="build", help="Genome build")
parser.add_argument("--project-id", dest="project_id", help="Google Project ID")

args = parser.parse_args()
vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores
mem = int(np.floor(float(args.mem)))
build = args.build
gcp_project = args.project_id

hl.init(
    min_block_size=128,
    spark_conf={
        "spark.executor.cores": cores,
        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
        "spark.driver.cores": cores,
        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
    },
    tmp_dir="tmp",
    local_tmpdir="tmp",
)


header = hl.get_vcf_metadata(vcf_file)

mt = hl.import_vcf(
    vcf_file,
    force_bgz=True,
    array_elements_required=False,
    call_fields=[],
    reference_genome=build,
)

mt = hl.vep(mt, config="vep_config.json", csq=True, tolerate_parse_error=True)

mt = mt.annotate_rows(info=mt.info.annotate(CSQ=mt.vep))

header["info"]["CSQ"] = {
    "Description": hl.eval(mt.vep_csq_header),
    "Number": ".",
    "Type": "String",
}

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)
