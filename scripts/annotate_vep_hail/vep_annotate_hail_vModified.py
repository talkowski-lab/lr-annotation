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

print(f"Reading VCF: {vcf_file}")
header = hl.get_vcf_metadata(vcf_file)
mt = hl.import_vcf(
    vcf_file,
    force_bgz=True,
    array_elements_required=False,
    call_fields=[],
    reference_genome=build,
)

print("Running VEP...")
mt = hl.vep(mt, config="vep_config.json", csq=False, tolerate_parse_error=True)

print("Selecting most severe consequence...")
most_severe_csq = mt.vep.most_severe_consequence

print("Flattening structure...")


def csq_struct_to_string(csq_struct):
    csq_fields = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "BIOTYPE",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "ALLELE_NUM",
        "DISTANCE",
        "STRAND",
        "FLAGS",
        "VARIANT_CLASS",
        "SYMBOL_SOURCE",
        "HGNC_ID",
        "CANONICAL",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "TSL",
        "APPRIS",
        "CCDS",
        "ENSP",
        "UNIPROT_ISOFORM",
        "SOURCE",
        "SIFT",
        "PolyPhen",
        "DOMAINS",
        "miRNA",
        "HGVS_OFFSET",
        "PUBMED",
        "MOTIF_NAME",
        "MOTIF_POS",
        "HIGH_INF_POS",
        "MOTIF_SCORE_CHANGE",
        "TRANSCRIPTION_FACTORS",
        "LoF",
        "LoF_filter",
        "LoF_flags",
        "LoF_info",
    ]

    fields = {
        "Allele": csq_struct.variant_allele,
        "Consequence": hl.delimit(csq_struct.consequence_terms, "&"),
        "IMPACT": csq_struct.impact,
        "SYMBOL": csq_struct.gene_symbol,
        "Gene": csq_struct.gene_id,
        "Feature_type": csq_struct.biotype,
        "Feature": csq_struct.transcript_id,
        "BIOTYPE": csq_struct.biotype,
        "CANONICAL": hl.if_else(csq_struct.canonical == 1, "YES", ""),
        "HGVSc": csq_struct.hgvsc,
        "HGVSp": csq_struct.hgvsp,
        "VARIANT_CLASS": mt.vep.variant_class,
        "SIFT": csq_struct.sift_prediction,
        "PolyPhen": csq_struct.polyphen_prediction,
    }

    return hl.delimit(
        [
            hl.or_else(hl.str(fields[f]) if f in fields else hl.null(hl.tstr), "")
            for f in csq_fields
        ],
        "|",
    )


print("Filtering for the canonical transcript consequence...")
transcript_csqs = hl.or_else(
    mt.vep.transcript_consequences,
    hl.empty_array(mt.vep.transcript_consequences.dtype.element_type),
)
canonical_csqs = transcript_csqs.filter(lambda csq: csq.canonical == 1)
csq = hl.if_else(
    hl.len(canonical_csqs) > 0,
    canonical_csqs[0],
    hl.null(canonical_csqs.dtype.element_type),
)

print("Formatting VEP string...")
vep_string = hl.if_else(
    hl.is_missing(csq),
    hl.delimit(
        ["", mt.vep.most_severe_consequence, "MODIFIER", "", "", "Intergenic", ""], "|"
    ),
    csq_struct_to_string(csq),
)
mt = mt.annotate_rows(info=mt.info.annotate(vep=vep_string))
mt = mt.drop("vep", "vep_proc_id")

print("Updating VCF header...")
gnomad_csq_header = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
    "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
    "Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|"
    "SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|"
    "APPRIS|CCDS|ENSP|UNIPROT_ISOFORM|SOURCE|SIFT|PolyPhen|DOMAINS|miRNA|"
    "HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|"
    "MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info"
)
header["info"]["vep"] = {
    "Description": (
        f"Consequence annotations from Ensembl VEP. Format: {gnomad_csq_header}"
    ),
    "Number": ".",
    "Type": "String",
}

print(f"Exporting annotated VCF to: {vep_annotated_vcf_name}")
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header)

print("VEP annotation complete.")
