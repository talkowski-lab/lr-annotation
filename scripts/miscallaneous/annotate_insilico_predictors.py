import argparse
import hail as hl

parser = argparse.ArgumentParser(description="Annotate VCF with GnomAD in silico predictors using Hail HTs.")

parser.add_argument('--build', required=True, help='Reference genome build (e.g., GRCh37 or GRCh38).')

parser.add_argument('--cadd_ht', required=True, help='CADD Hail Table path.')
parser.add_argument('--pangolin_ht', required=True, help='Pangolin Hail Table path.')
parser.add_argument('--phylop_ht', required=True, help='PhyloP Hail Table path.')
parser.add_argument('--revel_ht', required=True, help='REVEL Hail Table path.')
parser.add_argument('--spliceai_ht', required=True, help='SpliceAI Hail Table path.')

parser.add_argument('--vcf', required=True, help='Input VCF path.')
parser.add_argument('--output_vcf', required=True, help='Output annotated VCF path.')

args = parser.parse_args()

# Reference genome
build = args.build

# Annotation HT paths
cadd_ht_uri = args.cadd_ht
pangolin_ht_uri = args.pangolin_ht
phylop_ht_uri = args.phylop_ht
revel_ht_uri = args.revel_ht
spliceai_ht_uri = args.spliceai_ht

# VCF paths
vcf_uri = args.vcf
output_vcf_uri = args.output_vcf

hl.init(default_reference=build,
        min_block_size=128, 
        local=f"local[*]", 
        tmp_dir="tmp", local_tmpdir="tmp"
)

# Load all predictor HTs
cadd_ht = hl.read_table(cadd_ht_uri)
pangolin_ht = hl.read_table(pangolin_ht_uri)
phylop_ht = hl.read_table(phylop_ht_uri)
revel_ht = hl.read_table(revel_ht_uri)
spliceai_ht = hl.read_table(spliceai_ht_uri)

# Extract predictor versions from HT metadata
cadd_version = hl.eval(cadd_ht.cadd_version)
pangolin_version = hl.eval(pangolin_ht.pangolin_version)
phylop_version = hl.eval(phylop_ht.phylop_version)
revel_version = hl.eval(revel_ht.revel_version)
spliceai_version = hl.eval(spliceai_ht.spliceai_version)

# Import VCF
mt = hl.import_vcf(vcf_uri, force_bgz=True, array_elements_required=False, reference_genome=build)

# Get VCF header
header = hl.get_vcf_metadata(vcf_uri)

# Annotate rows with predictor scores
mt = mt.annotate_rows(
    info=mt.info.annotate(
        CADD_raw_score=cadd_ht[mt.locus, mt.alleles].cadd.raw_score,
        CADD_PHRED_score=cadd_ht[mt.locus, mt.alleles].cadd.phred,
        pangolin_largest_ds=pangolin_ht[mt.locus, mt.alleles].pangolin_largest_ds,
        revel_max=revel_ht[mt.locus, mt.alleles].revel_max,
        phylop=phylop_ht[mt.locus].phylop,
        spliceai_ds_max=spliceai_ht[mt.locus, mt.alleles].spliceai_ds_max
    )
)

# Predictor info and versions
predictors = {
    'CADD_raw_score': f"CADD {cadd_version}",
    'CADD_PHRED_score': f"CADD {cadd_version}",
    'pangolin_largest_ds': f"Pangolin versions {' and '.join(pangolin_version)}",
    'revel_max': f"REVEL {revel_version}",
    'phylop': f"PhyloP {phylop_version}",
    'spliceai_ds_max': f"SpliceAI {spliceai_version}"
}

# Update VCF header in a loop
for key, desc in predictors.items():
    header['info'][key] = {'Description': f"{key.replace('_', ' ').title()} from {desc}.", 'Number': '.', 'Type': 'Float'}

hl.export_vcf(dataset=mt, output=output_vcf_uri, metadata=header, tabix=True)
