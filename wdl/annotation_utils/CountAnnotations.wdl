version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "IntegrateVcfs.wdl" as IntegrateVcfs

workflow CountAnnotations {
	input {
		Array[File] vcfs
		Array[File] vcf_idxs
		String prefix

		Boolean create_per_sample = false
		Boolean create_per_allele = false
		Boolean create_list = false
		
		Boolean split_by_region = false
		Boolean create_variant_attributes = false
		String subset_vcf_string = ""
		Int max_length = -1
		Int min_length = -1
		
		Int? records_per_shard
		Int? shard_bin_size
		File? ref_fai

		String utils_docker

		RuntimeAttr? runtime_attr_subset
		RuntimeAttr? runtime_attr_shard
		RuntimeAttr? runtime_attr_create_shards
		RuntimeAttr? runtime_attr_region_subset
		RuntimeAttr? runtime_attr_count
		RuntimeAttr? runtime_attr_merge
	}

	scatter (i in range(length(vcfs))) {
		if (defined(shard_bin_size)) {
			call Helpers.CreateShardsFromVcfIndex {
				input:
					vcf_idx = vcf_idxs[i],
					ref_fai = select_first([ref_fai]),
					shard_bin_size = select_first([shard_bin_size]),
					prefix = "~{prefix}.input_~{i}.shards",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_create_shards
			}

			scatter (k in range(length(CreateShardsFromVcfIndex.shard_regions))) {
				call Helpers.SubsetVcfToRegionStreaming {
					input:
						vcf = vcfs[i],
						vcf_idx = vcf_idxs[i],
						region = CreateShardsFromVcfIndex.shard_regions[k],
						prefix = "~{prefix}.input_~{i}.region_~{k}",
						docker = utils_docker,
						runtime_attr_override = runtime_attr_region_subset
				}
			}
		}

		if (!defined(shard_bin_size) && defined(records_per_shard)) {
			call Helpers.ShardVcfByRecords {
				input:
					vcf = vcfs[i],
					vcf_idx = vcf_idxs[i],
					records_per_shard = select_first([records_per_shard]),
					prefix = "~{prefix}.input_~{i}",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_shard
			}
		}

		Array[File] shard_vcfs = select_first([SubsetVcfToRegionStreaming.subset_vcf, ShardVcfByRecords.shards, [vcfs[i]]])
		Array[File] shard_vcf_idxs = select_first([SubsetVcfToRegionStreaming.subset_vcf_idx, ShardVcfByRecords.shard_idxs, [vcf_idxs[i]]])

		scatter (j in range(length(shard_vcfs))) {
			if (create_variant_attributes) {
				call IntegrateVcfs.AnnotateVariantAttributes as AnnotateShardVariantAttributes {
					input:
						vcf = shard_vcfs[j],
						vcf_idx = shard_vcf_idxs[j],
						prefix = "~{prefix}.input_~{i}.shard_~{j}.annotated",
						docker = utils_docker,
						runtime_attr_override = runtime_attr_count
				}
			}
			
			call Helpers.SubsetVcfByArgs {
				input:
					vcf = select_first([AnnotateShardVariantAttributes.annotated_vcf, shard_vcfs[j]]),
					vcf_idx = select_first([AnnotateShardVariantAttributes.annotated_vcf_idx, shard_vcf_idxs[j]]),
					extra_args = subset_vcf_string,
					prefix = "~{prefix}.input_~{i}.shard_~{j}.subset",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset
			}

			call CountAnnotationShard {
				input:
					vcf = SubsetVcfByArgs.subset_vcf,
					vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
					create_per_sample = create_per_sample,
					create_per_allele = create_per_allele,
					create_list = create_list,
					split_by_region = split_by_region,
					max_length = max_length,
					min_length = min_length,
					prefix = "~{prefix}.input_~{i}.shard_~{j}",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_count
			}
		}
	}

	Array[File] site_count_tables = flatten(CountAnnotationShard.site_counts_tsv)
	Array[File] sample_count_tables = flatten(CountAnnotationShard.sample_counts_tsv)
	Array[File] allele_count_tables = flatten(CountAnnotationShard.allele_counts_tsv)
	Array[File] list_tables = flatten(CountAnnotationShard.list_tsv)
	Array[File] sample_count_files = flatten(CountAnnotationShard.sample_count_file)
	Array[File] functional_site_count_tables = flatten(CountAnnotationShard.gene_counts_tsv)
	Array[File] functional_sample_count_tables = flatten(CountAnnotationShard.gene_sample_counts_tsv)
	Array[File] functional_allele_count_tables = flatten(CountAnnotationShard.gene_allele_counts_tsv)

	call MergeAnnotationCountTables as MergeSiteCounts {
		input:
			count_tsvs = site_count_tables,
			sample_count_files = sample_count_files,
			normalization_mode = "sites",
			split_by_region = split_by_region,
			prefix = "~{prefix}.counts_sites",
			docker = utils_docker,
			runtime_attr_override = runtime_attr_merge
	}

	call MergeGeneCountTables as MergeGeneCounts {
		input:
			count_tsvs = functional_site_count_tables,
			prefix = "~{prefix}.counts_functional",
			docker = utils_docker,
			runtime_attr_override = runtime_attr_merge
	}

	if (create_per_sample) {
		call MergeGeneCountTables as MergeGeneSampleCounts {
			input:
				count_tsvs = functional_sample_count_tables,
				prefix = "~{prefix}.counts_functional_samples",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}

		call MergeAnnotationCountTables as MergeSampleCounts {
			input:
				count_tsvs = sample_count_tables,
				sample_count_files = sample_count_files,
				normalization_mode = "samples",
				split_by_region = split_by_region,
				prefix = "~{prefix}.counts_samples",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}
	}

	if (create_per_allele) {
		call MergeGeneCountTables as MergeGeneAlleleCounts {
			input:
				count_tsvs = functional_allele_count_tables,
				prefix = "~{prefix}.counts_functional_alleles",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}

		call MergeAnnotationCountTables as MergeAlleleCounts {
			input:
				count_tsvs = allele_count_tables,
				sample_count_files = sample_count_files,
				normalization_mode = "alleles",
				split_by_region = split_by_region,
				prefix = "~{prefix}.counts_alleles",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}
	}

	if (create_list) {
		call MergeAnnotationListTables {
			input:
				list_tsvs = list_tables,
				prefix = "~{prefix}.counts_list",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}
	}

	output {
		File annotation_counts_sites_tsv = MergeSiteCounts.merged_counts_tsv
		File annotation_counts_functional_tsv = MergeGeneCounts.merged_counts_tsv
		File? annotation_counts_functional_samples_tsv = MergeGeneSampleCounts.merged_counts_tsv
		File? annotation_counts_functional_alleles_tsv = MergeGeneAlleleCounts.merged_counts_tsv
		File? annotation_counts_list_tsv = MergeAnnotationListTables.merged_list_tsv
		File? annotation_counts_samples_tsv = MergeSampleCounts.merged_counts_tsv
		File? annotation_counts_alleles_tsv = MergeAlleleCounts.merged_counts_tsv
	}
}

task CountAnnotationShard {
	input {
		File vcf
		File vcf_idx
		Boolean create_per_sample
		Boolean create_per_allele
		Boolean create_list
		Boolean split_by_region
		Int max_length
		Int min_length
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<'PYCODE'
import csv
import re
from collections import defaultdict
import pysam

VCF_PATH = "~{vcf}"
SITE_OUTPUT = "~{prefix}.sites.raw.tsv"
SAMPLE_OUTPUT = "~{prefix}.samples.raw.tsv"
ALLELE_OUTPUT = "~{prefix}.alleles.raw.tsv"
LIST_OUTPUT = "~{prefix}.list.raw.tsv"
GENE_OUTPUT = "~{prefix}.genes.raw.tsv"
GENE_SAMPLE_OUTPUT = "~{prefix}.genes.samples.raw.tsv"
GENE_ALLELE_OUTPUT = "~{prefix}.genes.alleles.raw.tsv"
SAMPLE_COUNT_OUTPUT = "~{prefix}.sample_count.txt"
CREATE_PER_SAMPLE = "~{create_per_sample}".lower() == "true"
CREATE_PER_ALLELE = "~{create_per_allele}".lower() == "true"
CREATE_LIST = "~{create_list}".lower() == "true"
SPLIT_BY_REGION = "~{split_by_region}".lower() == "true"
MAX_LENGTH = ~{max_length}
MIN_LENGTH = ~{min_length}
MAX_LENGTH = MAX_LENGTH if MAX_LENGTH > 0 else None
MIN_LENGTH = MIN_LENGTH if MIN_LENGTH > 0 else None

COLUMN_BUCKETS = [
	"SNV",
	"INS 1-49bp",
	"DEL 1-49bp",
	"INS 50-499bp",
	"DEL 50-499bp",
	"INS >499bp",
	"DEL >499bp",
	"TRV",
	"Other",
]

INTERNAL_TOTAL_LABEL = "Total"
DISPLAY_ALL_LABEL = "All"
REGION_ORDER = ["US", "RM", "SD", "SR"]

ROW_ORDER = [
	("", INTERNAL_TOTAL_LABEL),
	("", "ME"),
	("ME", "ALU"),
	("ME", "LINE"),
	("ME", "SVA"),
	("", "Duplication"),
	("Duplication", "Tandem"),
	("Duplication", "Interspersed"),
	("Duplication", "Complex"),
	("Duplication", "NUMT"),
	("Duplication", "TR Parsed"),
	("", "Consequence"),
	("Consequence", "LoF (SVAnnotate)"),
	("Consequence", "LoF (VEP)"),
	("Consequence", "Missense"),
	("Consequence", "Coding"),
	("Consequence", "Intronic"),
	("Consequence", "Intergenic"),
	("Consequence", "Other"),
	("", "SVAnnotate"),
	("SVAnnotate", "Copy Gain"),
	("SVAnnotate", "IED"),
	("SVAnnotate", "PED"),
	("SVAnnotate", "TSSD"),
	("SVAnnotate", "DP"),
	("SVAnnotate", "UTR"),
	("", "Matched"),
	("Matched", "dbSNP Matched"),
	("Matched", "dbSNP Missing"),
	("Matched", "dbSNP Missing + LoF (SVAnnotate)"),
	("Matched", "dbSNP Missing + LoF (VEP)"),
	("Matched", "gnomAD Matched"),
	("Matched", "gnomAD Missing"),
	("Matched", "gnomAD Missing + LoF (SVAnnotate)"),
	("Matched", "gnomAD Missing + LoF (VEP)"),
	("Matched", "dbSNP/gnomAD Matched"),
	("Matched", "dbSNP/gnomAD Missing"),
	("Matched", "dbSNP/gnomAD Missing + LoF (SVAnnotate)"),
	("Matched", "dbSNP/gnomAD Missing + LoF (VEP)"),
]

PARENT_ROWS = {"ME", "Duplication", "Consequence", "SVAnnotate", "Matched"}

CONSEQUENCE_PRIORITY = [
	("LoF (VEP)", {
		"transcript_ablation", "stop_gained", "frameshift_variant",
		"splice_donor_variant", "splice_acceptor_variant", "stop_lost",
		"start_lost", "transcript_amplification", "feature_elongation",
		"feature_truncation",
	}),
	("Missense", {
		"missense_variant", "inframe_insertion", "inframe_deletion",
		"protein_altering_variant",
	}),
	("Coding", {
		"synonymous_variant", "stop_retained_variant", "start_retained_variant",
		"incomplete_terminal_codon_variant", "coding_sequence_variant",
	}),
	("Intronic", {
		"intron_variant", "splice_region_variant", "splice_donor_region_variant",
		"splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant",
		"NMD_transcript_variant", "non_coding_transcript_exon_variant",
		"non_coding_transcript_variant",
	}),
	("Intergenic", {
		"intergenic_variant", "upstream_gene_variant", "downstream_gene_variant",
		"regulatory_region_variant", "regulatory_region_ablation",
		"regulatory_region_amplification", "TF_binding_site_variant",
		"TFBS_ablation", "TFBS_amplification", "3_prime_UTR_variant",
		"5_prime_UTR_variant",
	}),
]

gene_counts = defaultdict(lambda: defaultdict(int))
gene_sample_counts = defaultdict(lambda: defaultdict(int))
gene_allele_counts = defaultdict(lambda: defaultdict(int))
PREDICTED_FIELDS = [
	"PREDICTED_LOF",
	"PREDICTED_COPY_GAIN",
	"PREDICTED_INTRAGENIC_EXON_DUP",
	"PREDICTED_PARTIAL_EXON_DUP",
	"PREDICTED_TSS_DUP",
	"PREDICTED_DUP_PARTIAL",
	"PREDICTED_INV_SPAN",
	"PREDICTED_UTR",
	"PREDICTED_INTRONIC",
	"PREDICTED_PROMOTER"
]

def init_table():
	return defaultdict(lambda: {column: 0.0 for column in COLUMN_BUCKETS})

def first_value(value):
	if isinstance(value, (list, tuple)):
		return value[0] if value else None
	return value

def get_info_value(record, key):
	try:
		return record.info.get(key)
	except ValueError:
		return None

def get_string_info(record, key):
	value = first_value(get_info_value(record, key))
	return "" if value is None else str(value)

def get_int_info(record, key):
	value = first_value(get_info_value(record, key))
	if value is None or value == ".":
		return None
	return abs(int(value))

def get_float_info(record, key):
	value = first_value(get_info_value(record, key))
	if value is None or value == ".":
		return None
	return float(value)

def has_info(record, key):
	return key in record.info

def get_info_gene_names(record, key):
	value = get_info_value(record, key)
	if isinstance(value, str):
		values = [value]
	elif isinstance(value, (list, tuple)):
		values = value
	elif value is None:
		values = []
	else:
		values = [str(value)]

	genes = set()
	for item in values:
		if item is None:
			continue
		for gene_name in str(item).split(","):
			gene_name = gene_name.strip()
			if gene_name and gene_name != ".":
				genes.add(gene_name)
	return genes

def get_vep_field_indices(header):
	if "vep" not in header.info:
		return {}
	description = header.info["vep"].description or ""
	match = re.search(r"Format:\s*([^\"]+)", description)
	if match is None:
		return {}
	fields = [field.strip() for field in match.group(1).strip().rstrip(".").split("|")]
	return {field: idx for idx, field in enumerate(fields)}

def extract_all_consequences(record, vep_field_indices):
	consequence_idx = vep_field_indices.get("Consequence")
	if consequence_idx is None or "vep" not in record.info:
		return set()
	annotations = get_info_value(record, "vep")
	if isinstance(annotations, str):
		annotations = [annotations]
	consequences = set()
	for annotation in annotations:
		fields = annotation.split("|")
		if consequence_idx >= len(fields):
			continue
		for consequence in fields[consequence_idx].split("&"):
			consequence = consequence.strip()
			if consequence:
				consequences.add(consequence)
	return consequences

def determine_consequence_labels(record, consequences):
	labels = set()
	if has_info(record, "PREDICTED_LOF"):
		labels.add("LoF (SVAnnotate)")
	for label, consequence_terms in CONSEQUENCE_PRIORITY:
		if consequence_terms & consequences:
			labels.add(label)
	if not labels:
		labels.add("Other")
	return labels

def determine_column(record):
	allele_type = get_string_info(record, "allele_type").lower()
	allele_length = get_int_info(record, "allele_length")
	variant_id = (record.id or "").upper()
	if allele_type == "snv": return "SNV"
	if (allele_type == "ins" or "INS" in variant_id) and allele_length < 50: return "INS 1-49bp"
	if (allele_type == "del" or "DEL" in variant_id) and allele_length < 50: return "DEL 1-49bp"
	if (allele_type == "ins" or "INS" in variant_id) and allele_length < 500: return "INS 50-499bp"
	if (allele_type == "del" or "DEL" in variant_id) and allele_length < 500: return "DEL 50-499bp"
	if (allele_type == "ins" or "INS" in variant_id) and allele_length >= 500: return "INS >499bp"
	if (allele_type == "del" or "DEL" in variant_id) and allele_length >= 500: return "DEL >499bp"
	if allele_type == "trv": return "TRV"
	return "Other"

def determine_row_weights(record, vep_field_indices):
	allele_type = get_string_info(record, "allele_type").lower()
	consequences = extract_all_consequences(record, vep_field_indices)
	consequence_labels = determine_consequence_labels(record, consequences)

	row_weights = {row: 0 for row in ROW_ORDER}
	row_weights[("", INTERNAL_TOTAL_LABEL)] = 1

	is_alu = "alu" in allele_type
	is_line = "line" in allele_type
	is_sva = "sva" in allele_type
	if is_alu or is_line or is_sva: row_weights[("", "ME")] = 1
	if is_alu: row_weights[("ME", "ALU")] = 1
	if is_line: row_weights[("ME", "LINE")] = 1
	if is_sva: row_weights[("ME", "SVA")] = 1

	is_dup_interspersed = "dup_interspersed" in allele_type
	is_dup_complex = "complex_dup" in allele_type
	is_dup_tandem = "dup" in allele_type and not is_dup_interspersed and not is_dup_complex
	is_numt = "numt" in allele_type
	is_tr_parsed = has_info(record, "TR_PARSED")
	
	if is_dup_tandem or is_dup_interspersed or is_dup_complex or is_numt or is_tr_parsed:
		row_weights[("", "Duplication")] = 1
	if is_dup_tandem: row_weights[("Duplication", "Tandem")] = 1
	if is_dup_interspersed: row_weights[("Duplication", "Interspersed")] = 1
	if is_dup_complex: row_weights[("Duplication", "Complex")] = 1
	if is_numt: row_weights[("Duplication", "NUMT")] = 1
	if is_tr_parsed: row_weights[("Duplication", "TR Parsed")] = 1

	row_weights[("", "Consequence")] = 1
	for consequence_label in consequence_labels:
		row_weights[("Consequence", consequence_label)] = 1

	is_copy_gain = has_info(record, "PREDICTED_COPY_GAIN")
	is_ied = has_info(record, "PREDICTED_INTRAGENIC_EXON_DUP")
	is_ped = has_info(record, "PREDICTED_PARTIAL_EXON_DUP")
	is_tssd = has_info(record, "PREDICTED_TSS_DUP")
	is_dp = has_info(record, "PREDICTED_DUP_PARTIAL")
	is_utr = has_info(record, "PREDICTED_UTR")
	is_dbsnp = has_info(record, "dbSNP_ID") or has_info(record, "dbGaP_ID")
	is_gnomad_matched = has_info(record, "gnomAD_V4_match_ID")
	has_svannotate = is_copy_gain or is_ied or is_ped or is_tssd or is_dp or is_utr

	if has_svannotate: row_weights[("", "SVAnnotate")] = 1
	if is_copy_gain: row_weights[("SVAnnotate", "Copy Gain")] = 1
	if is_ied: row_weights[("SVAnnotate", "IED")] = 1
	if is_ped: row_weights[("SVAnnotate", "PED")] = 1
	if is_tssd: row_weights[("SVAnnotate", "TSSD")] = 1
	if is_dp: row_weights[("SVAnnotate", "DP")] = 1
	if is_utr: row_weights[("SVAnnotate", "UTR")] = 1

	row_weights[("", "Matched")] = 1
	if is_dbsnp:
		row_weights[("Matched", "dbSNP Matched")] = 1
	else:
		row_weights[("Matched", "dbSNP Missing")] = 1

	if is_gnomad_matched:
		row_weights[("Matched", "gnomAD Matched")] = 1
	else:
		row_weights[("Matched", "gnomAD Missing")] = 1

	if is_dbsnp or is_gnomad_matched:
		row_weights[("Matched", "dbSNP/gnomAD Matched")] = 1
	else:
		row_weights[("Matched", "dbSNP/gnomAD Missing")] = 1

	if "LoF (SVAnnotate)" in consequence_labels:
		if not is_dbsnp:
			row_weights[("Matched", "dbSNP Missing + LoF (SVAnnotate)")] = 1
		if not is_gnomad_matched:
			row_weights[("Matched", "gnomAD Missing + LoF (SVAnnotate)")] = 1
		if not is_dbsnp and not is_gnomad_matched:
			row_weights[("Matched", "dbSNP/gnomAD Missing + LoF (SVAnnotate)")] = 1
	if "LoF (VEP)" in consequence_labels:
		if not is_dbsnp:
			row_weights[("Matched", "dbSNP Missing + LoF (VEP)")] = 1
		if not is_gnomad_matched:
			row_weights[("Matched", "gnomAD Missing + LoF (VEP)")] = 1
		if not is_dbsnp and not is_gnomad_matched:
			row_weights[("Matched", "dbSNP/gnomAD Missing + LoF (VEP)")] = 1

	return {row_key: weight for row_key, weight in row_weights.items() if weight > 0}

def get_genotype_weights(record):
	carrier_count = 0
	alt_allele_count = 0
	for sample in record.samples.values():
		genotype = sample.get("GT")
		if genotype is None: continue
		alt_alleles = sum(1 for allele in genotype if allele is not None and allele > 0)
		alt_allele_count += alt_alleles
		if alt_alleles > 0: carrier_count += 1
	return carrier_count, alt_allele_count

def format_list_column(row_key):
	annotation_group, annotation = row_key
	if annotation_group: return annotation
	if annotation in PARENT_ROWS: return f"{annotation} - {INTERNAL_TOTAL_LABEL}"
	return annotation

def get_list_columns():
	default_labels = [format_list_column(row_key) for row_key in ROW_ORDER]
	label_counts = {}
	for label in default_labels:
		label_counts[label] = label_counts.get(label, 0) + 1
	list_columns = []
	for row_key, default_label in zip(ROW_ORDER, default_labels):
		annotation_group, annotation = row_key
		if annotation_group and label_counts[default_label] > 1:
			list_columns.append(f"{annotation_group} - {annotation}")
		else:
			list_columns.append(default_label)
	return list_columns

def write_table(path, table_data, integer_output, split_by_region):
	with open(path, "w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		header = ["category", "sub_category", "tr_status"]
		if split_by_region:
			header.append("region")
		writer.writerow(header + COLUMN_BUCKETS)
		for key in sorted(table_data.keys()):
			if split_by_region:
				category, sub_category, tr_status, region = key
			else:
				category, sub_category, tr_status = key
			display_category = DISPLAY_ALL_LABEL if category == INTERNAL_TOTAL_LABEL else category
			display_sub_category = DISPLAY_ALL_LABEL if sub_category == INTERNAL_TOTAL_LABEL else sub_category
			display_tr_status = DISPLAY_ALL_LABEL if tr_status == INTERNAL_TOTAL_LABEL else tr_status
			row_prefix = [display_category, display_sub_category, display_tr_status]
			if split_by_region:
				display_region = DISPLAY_ALL_LABEL if region == INTERNAL_TOTAL_LABEL else region
				row_prefix.append(display_region)
			values = []
			for column in COLUMN_BUCKETS:
				value = table_data[key][column]
				if integer_output: values.append(str(int(value)))
				else: values.append(str(value))
			writer.writerow(row_prefix + values)

def write_gene_counts(path, gene_count_data):
	with open(path, "w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		writer.writerow(["gene_name", "category", "count"])
		for gene in gene_count_data:
			for category, count in gene_count_data[gene].items():
				writer.writerow([gene, category, str(count)])

site_table = init_table()
sample_table = init_table()
allele_table = init_table()

vcf_in = pysam.VariantFile(VCF_PATH)
vep_field_indices = get_vep_field_indices(vcf_in.header)
sample_count = len(vcf_in.header.samples)

with open(SAMPLE_COUNT_OUTPUT, "w") as handle:
	handle.write(f"{sample_count}\n")

with open(LIST_OUTPUT, "w", newline="") as handle:
	writer = csv.writer(handle, delimiter="\t")
	writer.writerow(["variant_id", "classification"] + get_list_columns())

	for record in vcf_in:
		allele_length = get_int_info(record, "allele_length")
		if MAX_LENGTH is not None and (allele_length is None or allele_length > MAX_LENGTH):
			continue
		if MIN_LENGTH is not None and (allele_length is None or allele_length < MIN_LENGTH):
			continue

		column = determine_column(record)
		row_weights = determine_row_weights(record, vep_field_indices)
		
		tr_status = "TR" if has_info(record, "TR_ENVELOPED") else "Not TR"
		regions = [INTERNAL_TOTAL_LABEL]
		if SPLIT_BY_REGION:
			region = get_string_info(record, "REGION")
			if region in REGION_ORDER:
				regions.append(region)

		if CREATE_PER_SAMPLE or CREATE_PER_ALLELE:
			carrier_count, alt_allele_count = get_genotype_weights(record)
		else:
			carrier_count, alt_allele_count = 0, 0

		af = get_float_info(record, "AF")
		ac = get_int_info(record, "AC")
		
		is_large = allele_length is not None and allele_length > 50000
		
		buckets = ["all"]
		if is_large:
			buckets.append("all_large")
		if ac == 1:
			buckets.append("singleton")
			if is_large:
				buckets.append("singleton_large")
		if af is not None and af < 0.001:
			buckets.append("ultra_rare")
			if is_large:
				buckets.append("ultra_rare_large")
		if af is not None and af < 0.01:
			buckets.append("rare")
			if is_large:
				buckets.append("rare_large")
		if af is not None and af > 0.05:
			buckets.append("common")
			if is_large:
				buckets.append("common_large")
		
		for field in PREDICTED_FIELDS:
			if has_info(record, field):
				for gene_name in get_info_gene_names(record, field):
					for bucket in buckets:
						category = f"{field}.{bucket}"
						gene_counts[gene_name][category] += 1
						if CREATE_PER_SAMPLE:
							gene_sample_counts[gene_name][category] += carrier_count
						if CREATE_PER_ALLELE:
							gene_allele_counts[gene_name][category] += alt_allele_count

		for row_key, weight in row_weights.items():
			cat = row_key[0] if row_key[0] else row_key[1]
			ann = row_key[1] if row_key[0] else INTERNAL_TOTAL_LABEL
			for current_tr_status in [INTERNAL_TOTAL_LABEL, tr_status]:
				if SPLIT_BY_REGION:
					for current_region in regions:
						full_key = (cat, ann, current_tr_status, current_region)
						site_table[full_key][column] += weight
						if CREATE_PER_SAMPLE: sample_table[full_key][column] += carrier_count * weight
						if CREATE_PER_ALLELE: allele_table[full_key][column] += alt_allele_count * weight
				else:
					full_key = (cat, ann, current_tr_status)
					site_table[full_key][column] += weight
					if CREATE_PER_SAMPLE: sample_table[full_key][column] += carrier_count * weight
					if CREATE_PER_ALLELE: allele_table[full_key][column] += alt_allele_count * weight

		if CREATE_LIST:
			writer.writerow([
				str(record.id or "."),
				column,
			] + ["1" if row_key in row_weights else "0" for row_key in ROW_ORDER])

write_table(SITE_OUTPUT, site_table, integer_output=True, split_by_region=SPLIT_BY_REGION)
write_table(SAMPLE_OUTPUT, sample_table, integer_output=True, split_by_region=SPLIT_BY_REGION)
write_table(ALLELE_OUTPUT, allele_table, integer_output=True, split_by_region=SPLIT_BY_REGION)
write_gene_counts(GENE_OUTPUT, gene_counts)
write_gene_counts(GENE_SAMPLE_OUTPUT, gene_sample_counts)
write_gene_counts(GENE_ALLELE_OUTPUT, gene_allele_counts)

if not CREATE_LIST:
	with open(LIST_OUTPUT, "w", newline="") as handle:
		pass
PYCODE
	>>>

	output {
		File site_counts_tsv = "~{prefix}.sites.raw.tsv"
		File sample_counts_tsv = "~{prefix}.samples.raw.tsv"
		File allele_counts_tsv = "~{prefix}.alleles.raw.tsv"
		File list_tsv = "~{prefix}.list.raw.tsv"
		File gene_counts_tsv = "~{prefix}.genes.raw.tsv"
		File gene_sample_counts_tsv = "~{prefix}.genes.samples.raw.tsv"
		File gene_allele_counts_tsv = "~{prefix}.genes.alleles.raw.tsv"
		File sample_count_file = "~{prefix}.sample_count.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
		boot_disk_gb: 10,
		preemptible_tries: 2,
		max_retries: 0
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: 8 + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task MergeAnnotationCountTables {
	input {
		Array[File] count_tsvs
		Array[File] sample_count_files
		String normalization_mode
		Boolean split_by_region = false
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<'PYCODE'
import csv

COUNT_FILES = "~{sep=',' count_tsvs}".split(",")
SAMPLE_COUNT_FILES = [path for path in "~{sep=',' sample_count_files}".split(",") if path]
MODE = "~{normalization_mode}"
SPLIT_BY_REGION = "~{split_by_region}".lower() == "true"
OUTPUT = "~{prefix}.tsv"

INTERNAL_TOTAL_LABEL = "Total"
DISPLAY_ALL_LABEL = "All"
TR_ORDER = [INTERNAL_TOTAL_LABEL, "TR", "Not TR"]
REGION_ORDER = [INTERNAL_TOTAL_LABEL, "US", "RM", "SD", "SR"]

ROW_ORDER = [
	("", INTERNAL_TOTAL_LABEL),
	("", "ME"),
	("ME", "ALU"),
	("ME", "LINE"),
	("ME", "SVA"),
	("", "Duplication"),
	("Duplication", "Tandem"),
	("Duplication", "Interspersed"),
	("Duplication", "Complex"),
	("Duplication", "NUMT"),
	("Duplication", "TR Parsed"),
	("", "Consequence"),
	("Consequence", "LoF (SVAnnotate)"),
	("Consequence", "LoF (VEP)"),
	("Consequence", "Missense"),
	("Consequence", "Coding"),
	("Consequence", "Intronic"),
	("Consequence", "Intergenic"),
	("Consequence", "Other"),
	("", "SVAnnotate"),
	("SVAnnotate", "Copy Gain"),
	("SVAnnotate", "IED"),
	("SVAnnotate", "PED"),
	("SVAnnotate", "TSSD"),
	("SVAnnotate", "DP"),
	("SVAnnotate", "UTR"),
	("", "Matched"),
	("Matched", "dbSNP Matched"),
	("Matched", "dbSNP Missing"),
	("Matched", "dbSNP Missing + LoF (SVAnnotate)"),
	("Matched", "dbSNP Missing + LoF (VEP)"),
	("Matched", "gnomAD Matched"),
	("Matched", "gnomAD Missing"),
	("Matched", "gnomAD Missing + LoF (SVAnnotate)"),
	("Matched", "gnomAD Missing + LoF (VEP)"),
	("Matched", "dbSNP/gnomAD Matched"),
	("Matched", "dbSNP/gnomAD Missing"),
	("Matched", "dbSNP/gnomAD Missing + LoF (SVAnnotate)"),
	("Matched", "dbSNP/gnomAD Missing + LoF (VEP)"),
]

def get_cat_ann(row_key):
	group, ann = row_key
	if not group:
		return ann, INTERNAL_TOTAL_LABEL
	return group, ann

def normalize_input_label(value):
	return INTERNAL_TOTAL_LABEL if value == DISPLAY_ALL_LABEL else value

header = None
counts = {}
group_column_count = 4 if SPLIT_BY_REGION else 3

for path in COUNT_FILES:
	with open(path, "r", newline="") as handle:
		reader = csv.reader(handle, delimiter="\t")
		current_header = next(reader)
		if header is None:
			header = current_header
		elif current_header != header:
			raise ValueError(f"Mismatched headers while merging count tables: {path}")

		for row in reader:
			key = tuple(normalize_input_label(value) for value in row[:group_column_count])
			if key not in counts:
				counts[key] = [0.0] * (len(header) - group_column_count)
			for i, val in enumerate(row[group_column_count:]):
				counts[key][i] += float(val)

if header is None:
	raise ValueError("No count tables were provided for merging")

denominator = 1.0
if MODE != "sites":
	sample_counts = set()
	for path in SAMPLE_COUNT_FILES:
		with open(path, "r") as handle:
			sample_counts.add(int(handle.read().strip()))

	if len(sample_counts) != 1:
		raise ValueError("All input VCFs must contain the same number of samples for normalized count outputs")

	sample_count = sample_counts.pop()
	if sample_count <= 0:
		raise ValueError("Sample-normalized outputs require at least one sample in every input VCF")

	denominator = float(sample_count)

value_count = len(header) - group_column_count

with open(OUTPUT, "w", newline="") as handle:
	writer = csv.writer(handle, delimiter="\t")
	writer.writerow(header)

	def get_row_prefix(c, a, t, region, suffix=""):
		display_category = DISPLAY_ALL_LABEL if c == INTERNAL_TOTAL_LABEL else c
		display_sub_category = DISPLAY_ALL_LABEL if a == INTERNAL_TOTAL_LABEL else a
		display_tr_status = DISPLAY_ALL_LABEL if t == INTERNAL_TOTAL_LABEL else t
		row_prefix = [display_category, f"{display_sub_category}{suffix}", display_tr_status]
		if SPLIT_BY_REGION:
			display_region = DISPLAY_ALL_LABEL if region == INTERNAL_TOTAL_LABEL else region
			row_prefix.append(display_region)
		return row_prefix

	def write_count_row(c, a, t, region, values):
		if MODE == "sites":
			formatted = [str(int(round(v))) for v in values]
		else:
			formatted = [f"{v / denominator:.2f}" for v in values]

		writer.writerow(get_row_prefix(c, a, t, region) + formatted)

	def write_percentage_row(c, a, t, region, values, denominator_values):
		formatted = []
		for value, denominator_value in zip(values, denominator_values):
			if denominator_value == 0:
				formatted.append("0.0000")
			else:
				formatted.append(f"{value / denominator_value:.4f}")

		writer.writerow(get_row_prefix(c, a, t, region, suffix=" (prop)") + formatted)

	for row_key in ROW_ORDER:
		cat, ann = get_cat_ann(row_key)
		for tr in TR_ORDER:
			if SPLIT_BY_REGION:
				for region in REGION_ORDER:
					values = counts.get((cat, ann, tr, region), [0.0] * value_count)
					denominator_values = counts.get((INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL, tr, region), [0.0] * value_count)
					write_count_row(cat, ann, tr, region, values)
					write_percentage_row(cat, ann, tr, region, values, denominator_values)
			else:
				values = counts.get((cat, ann, tr), [0.0] * value_count)
				denominator_values = counts.get((INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL, tr), [0.0] * value_count)
				write_count_row(cat, ann, tr, INTERNAL_TOTAL_LABEL, values)
				write_percentage_row(cat, ann, tr, INTERNAL_TOTAL_LABEL, values, denominator_values)
PYCODE
	>>>

	output {
		File merged_counts_tsv = "~{prefix}.tsv"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
		boot_disk_gb: 10,
		preemptible_tries: 2,
		max_retries: 0
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task MergeAnnotationListTables {
	input {
		Array[File] list_tsvs
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<'PYCODE'
import csv


LIST_FILES = "~{sep=',' list_tsvs}".split(",")
OUTPUT = "~{prefix}.tsv"


header = None
variant_order = []
variants = {}

for path in LIST_FILES:
	with open(path, "r", newline="") as handle:
		reader = csv.reader(handle, delimiter="\t")
		current_header = next(reader)
		if header is None:
			header = current_header
		elif current_header != header:
			raise ValueError(f"Mismatched headers while merging annotation list tables: {path}")

		for row in reader:
			variant_id = row[0]
			classification = row[1]
			memberships = [int(value) for value in row[2:]]

			if variant_id not in variants:
				variants[variant_id] = [classification, memberships]
				variant_order.append(variant_id)
				continue

			existing_classification, existing_memberships = variants[variant_id]
			if existing_classification != classification:
				raise ValueError(
					f"Variant {variant_id} has conflicting classifications while merging annotation list tables"
				)

			variants[variant_id][1] = [max(left, right) for left, right in zip(existing_memberships, memberships)]

if header is None:
	raise ValueError("No annotation list tables were provided for merging")

with open(OUTPUT, "w", newline="") as handle:
	writer = csv.writer(handle, delimiter="\t")
	writer.writerow(header)
	for variant_id in variant_order:
		classification, memberships = variants[variant_id]
		writer.writerow([variant_id, classification] + [str(value) for value in memberships])
PYCODE
	>>>

	output {
		File merged_list_tsv = "~{prefix}.tsv"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 4 * ceil(size(list_tsvs, "GB")) + 10,
		boot_disk_gb: 10,
		preemptible_tries: 2,
		max_retries: 0
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task MergeGeneCountTables {
	input {
		Array[File] count_tsvs
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<'PYCODE'
import csv
from collections import defaultdict

COUNT_FILES = "~{sep=',' count_tsvs}".split(",")
OUTPUT = "~{prefix}.tsv"

PREDICTED_FIELDS = [
	"PREDICTED_LOF",
	"PREDICTED_COPY_GAIN",
	"PREDICTED_INTRAGENIC_EXON_DUP",
	"PREDICTED_PARTIAL_EXON_DUP",
	"PREDICTED_TSS_DUP",
	"PREDICTED_DUP_PARTIAL",
	"PREDICTED_INV_SPAN",
	"PREDICTED_UTR",
	"PREDICTED_INTRONIC",
	"PREDICTED_PROMOTER"
]

BUCKETS = [
	"all",
	"common",
	"rare",
	"ultra_rare",
	"singleton",
	"all_large",
	"common_large",
	"rare_large",
	"ultra_rare_large",
	"singleton_large",
]

new_columns = []
for f in PREDICTED_FIELDS:
	for b in BUCKETS:
		new_columns.append(f"{f}.{b}")

gene_counts = defaultdict(lambda: {col: 0 for col in new_columns})

for path in COUNT_FILES:
	with open(path, "r", newline="") as handle:
		reader = csv.reader(handle, delimiter="\t")
		try:
			next(reader) # skip header
		except StopIteration:
			pass
		for row in reader:
			if len(row) == 3:
				gene, cat, count = row
				gene_counts[gene][cat] += int(count)

with open(OUTPUT, "w", newline="") as out_handle:
	writer = csv.writer(out_handle, delimiter="\t")
	
	# Write header
	writer.writerow(["gene_name"] + new_columns)
		
	# Write data
	for gene in sorted(gene_counts.keys()):
		counts = gene_counts[gene]
		row = [gene] + [str(counts.get(col, 0)) for col in new_columns]
		writer.writerow(row)

PYCODE
	>>>

	output {
		File merged_counts_tsv = "~{prefix}.tsv"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
		boot_disk_gb: 10,
		preemptible_tries: 2,
		max_retries: 0
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}
