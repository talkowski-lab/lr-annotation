version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CountAnnotations {
	input {
		Array[File] vcfs
		Array[File] vcfs_idx
		String prefix

		Boolean do_per_sample = false
		Boolean do_per_allele = false

		String? subset_vcf_string
		Int? records_per_shard

		String utils_docker

		RuntimeAttr? runtime_attr_subset
		RuntimeAttr? runtime_attr_shard
		RuntimeAttr? runtime_attr_count
		RuntimeAttr? runtime_attr_merge
	}

	scatter (i in range(length(vcfs))) {
		if (defined(subset_vcf_string)) {
			call SubsetVcfWithFlags {
				input:
					vcf = vcfs[i],
					vcf_idx = vcfs_idx[i],
					subset_flags = select_first([subset_vcf_string]),
					prefix = "~{prefix}.input_~{i}.subset",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset
			}
		}

		File vcf_to_process = select_first([SubsetVcfWithFlags.subset_vcf, vcfs[i]])
		File vcf_idx_to_process = select_first([SubsetVcfWithFlags.subset_vcf_idx, vcfs_idx[i]])

		if (defined(records_per_shard)) {
			call Helpers.ShardVcfByRecords {
				input:
					vcf = vcf_to_process,
					vcf_idx = vcf_idx_to_process,
					records_per_shard = select_first([records_per_shard]),
					prefix = "~{prefix}.input_~{i}",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_shard
			}
		}

		Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [vcf_to_process]])
		Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx_to_process]])

		scatter (j in range(length(shard_vcfs))) {
			call CountAnnotationShard {
				input:
					vcf = shard_vcfs[j],
					vcf_idx = shard_vcf_idxs[j],
					do_per_sample = do_per_sample,
					do_per_allele = do_per_allele,
					prefix = "~{prefix}.input_~{i}.shard_~{j}",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_count
			}
		}
	}

	Array[File] site_count_tables = flatten(CountAnnotationShard.site_counts_tsv)
	Array[File] sample_count_tables = flatten(CountAnnotationShard.sample_counts_tsv)
	Array[File] allele_count_tables = flatten(CountAnnotationShard.allele_counts_tsv)
	Array[File] sample_count_files = flatten(CountAnnotationShard.sample_count_file)

	call MergeAnnotationCountTables as MergeSiteCounts {
		input:
			count_tsvs = site_count_tables,
			sample_count_files = sample_count_files,
			normalization_mode = "sites",
			prefix = "~{prefix}.annotation_counts_sites",
			docker = utils_docker,
			runtime_attr_override = runtime_attr_merge
	}

	if (do_per_sample) {
		call MergeAnnotationCountTables as MergeSampleCounts {
			input:
				count_tsvs = sample_count_tables,
				sample_count_files = sample_count_files,
				normalization_mode = "samples",
				prefix = "~{prefix}.annotation_counts_samples",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}
	}

	if (do_per_allele) {
		call MergeAnnotationCountTables as MergeAlleleCounts {
			input:
				count_tsvs = allele_count_tables,
				sample_count_files = sample_count_files,
				normalization_mode = "alleles",
				prefix = "~{prefix}.annotation_counts_alleles",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_merge
		}
	}

	output {
		File annotation_counts_sites_tsv = MergeSiteCounts.merged_counts_tsv
		File? annotation_counts_samples_tsv = MergeSampleCounts.merged_counts_tsv
		File? annotation_counts_alleles_tsv = MergeAlleleCounts.merged_counts_tsv
	}
}

task SubsetVcfWithFlags {
	input {
		File vcf
		File vcf_idx
		String subset_flags
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		bcftools view \
			--no-version \
			~{subset_flags} \
			-O z \
			-o ~{prefix}.vcf.gz \
			~{vcf}

		tabix -p vcf ~{prefix}.vcf.gz
	>>>

	output {
		File subset_vcf = "~{prefix}.vcf.gz"
		File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task CountAnnotationShard {
	input {
		File vcf
		File vcf_idx
		Boolean do_per_sample
		Boolean do_per_allele
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<'PYCODE'
import csv
import re

import pysam


VCF_PATH = "~{vcf}"
SITE_OUTPUT = "~{prefix}.sites.raw.tsv"
SAMPLE_OUTPUT = "~{prefix}.samples.raw.tsv"
ALLELE_OUTPUT = "~{prefix}.alleles.raw.tsv"
SAMPLE_COUNT_OUTPUT = "~{prefix}.sample_count.txt"
DO_PER_SAMPLE = "~{do_per_sample}".lower() == "true"
DO_PER_ALLELE = "~{do_per_allele}".lower() == "true"

COLUMN_BUCKETS = [
	"SNV",
	"INS 1-49bp",
	"DEL 1-49bp",
	"INS >50bp",
	"DEL >50bp",
	"TRV",
	"Other",
]

ROW_ORDER = [
	("", "All"),
	("", "STR Overlap"),
	("", "ME"),
	("ME", "ALU"),
	("ME", "LINE"),
	("ME", "SVA"),
	("", "DUP"),
	("DUP", "DUP_TANDEM"),
	("DUP", "DUP_INTERSPERSED"),
	("DUP", "DUP_COMPLEX"),
	("", "NUMT"),
	("", "Intergenic"),
	("", "Intronic"),
	("", "LoF"),
	("", "Coding"),
	("", "SVAnnotate"),
	("SVAnnotate", "Copy Gain"),
	("SVAnnotate", "IED"),
	("SVAnnotate", "PED"),
	("SVAnnotate", "TSSD"),
	("SVAnnotate", "DP"),
	("SVAnnotate", "UTR"),
	("", "dbGaP"),
	("", "gnomAD Matched"),
]

LOF_CONSEQUENCES = {
	"transcript_ablation",
	"splice_acceptor_variant",
	"splice_donor_variant",
	"stop_gained",
	"frameshift_variant",
	"stop_lost",
	"start_lost",
	"transcript_amplification",
	"feature_elongation",
	"feature_truncation",
}

CODING_CONSEQUENCES = {
	"inframe_insertion",
	"inframe_deletion",
	"missense_variant",
	"protein_altering_variant",
}


def init_table():
	return {row: {column: 0.0 for column in COLUMN_BUCKETS} for row in ROW_ORDER}


def first_value(value):
	if isinstance(value, (list, tuple)):
		return value[0] if value else None
	return value


def get_string_info(record, key):
	value = first_value(record.info.get(key))
	return "" if value is None else str(value)


def get_int_info(record, key):
	value = first_value(record.info.get(key))
	if value is None or value == ".":
		return None
	return abs(int(value))


def has_info(record, key):
	return key in record.info


def get_vep_consequence_index(header):
	if "vep" not in header.info:
		return None
	description = header.info["vep"].description or ""
	match = re.search(r"Format:\s*([^\"]+)", description)
	if match is None:
		return None
	fields = [field.strip() for field in match.group(1).strip().rstrip(".").split("|")]
	for idx, field in enumerate(fields):
		if field == "Consequence":
			return idx
	return None


def extract_consequences(record, consequence_idx):
	if consequence_idx is None or "vep" not in record.info:
		return set()

	annotations = record.info.get("vep")
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


def determine_column(record):
	allele_type = get_string_info(record, "allele_type").lower()
	allele_length = get_int_info(record, "allele_length")
	variant_id = (record.id or "").upper()

	if allele_type == "snv":
		return "SNV"
	if allele_length is not None and "INS" in variant_id and allele_length < 50:
		return "INS 1-49bp"
	if allele_length is not None and "DEL" in variant_id and allele_length < 50:
		return "DEL 1-49bp"
	if allele_length is not None and "INS" in variant_id and allele_length >= 50:
		return "INS >50bp"
	if allele_length is not None and "DEL" in variant_id and allele_length >= 50:
		return "DEL >50bp"
	if allele_type == "trv":
		return "TRV"
	return "Other"


def determine_row_weights(record, consequence_idx):
	allele_type = get_string_info(record, "allele_type").lower()
	consequences = extract_consequences(record, consequence_idx)

	row_weights = {row: 0 for row in ROW_ORDER}
	row_weights[("", "All")] = 1

	if has_info(record, "TR_ENVELOPED"):
		row_weights[("", "STR Overlap")] = 1

	is_alu = "alu" in allele_type
	is_line = "line" in allele_type
	is_sva = "sva" in allele_type
	me_hits = int(is_alu) + int(is_line) + int(is_sva)
	if me_hits:
		row_weights[("", "ME")] = me_hits
	if is_alu:
		row_weights[("ME", "ALU")] = 1
	if is_line:
		row_weights[("ME", "LINE")] = 1
	if is_sva:
		row_weights[("ME", "SVA")] = 1

	is_dup_interspersed = "dup_interspersed" in allele_type
	is_dup_complex = "complex_dup" in allele_type
	is_dup_tandem = "dup" in allele_type and not is_dup_interspersed and not is_dup_complex
	dup_hits = int(is_dup_tandem) + int(is_dup_interspersed) + int(is_dup_complex)
	if dup_hits:
		row_weights[("", "DUP")] = dup_hits
	if is_dup_tandem:
		row_weights[("DUP", "DUP_TANDEM")] = 1
	if is_dup_interspersed:
		row_weights[("DUP", "DUP_INTERSPERSED")] = 1
	if is_dup_complex:
		row_weights[("DUP", "DUP_COMPLEX")] = 1

	if "numt" in allele_type:
		row_weights[("", "NUMT")] = 1

	is_intergenic = has_info(record, "PREDICTED_INTERGENIC") or "intergenic_variant" in consequences
	if is_intergenic:
		row_weights[("", "Intergenic")] = 1

	is_intronic = has_info(record, "PREDICTED_INTRONIC") or "intron_variant" in consequences
	if is_intronic:
		row_weights[("", "Intronic")] = 1

	is_lof = has_info(record, "PREDICTED_LOF") or bool(LOF_CONSEQUENCES & consequences)
	if is_lof:
		row_weights[("", "LoF")] = 1

	is_copy_gain = has_info(record, "PREDICTED_COPY_GAIN")
	is_ied = has_info(record, "PREDICTED_INTRAGENIC_EXON_DUP")
	is_ped = has_info(record, "PREDICTED_PARTIAL_EXON_DUP")
	is_tssd = has_info(record, "PREDICTED_TSS_DUP")
	is_dp = has_info(record, "PREDICTED_DUP_PARTIAL")
	is_utr = has_info(record, "PREDICTED_UTR")
	is_dbgap = has_info(record, "dbGaP_ID")
	is_gnomad_matched = has_info(record, "gnomAD_V4_match_ID")
	is_vep_coding = bool(CODING_CONSEQUENCES & consequences)
	svannotate_hits = int(is_copy_gain) + int(is_ied) + int(is_ped) + int(is_tssd) + int(is_dp) + int(is_utr)

	if svannotate_hits:
		row_weights[("", "SVAnnotate")] = svannotate_hits
	if is_copy_gain:
		row_weights[("SVAnnotate", "Copy Gain")] = 1
	if is_ied:
		row_weights[("SVAnnotate", "IED")] = 1
	if is_ped:
		row_weights[("SVAnnotate", "PED")] = 1
	if is_tssd:
		row_weights[("SVAnnotate", "TSSD")] = 1
	if is_dp:
		row_weights[("SVAnnotate", "DP")] = 1
	if is_utr:
		row_weights[("SVAnnotate", "UTR")] = 1

	coding_hits = int(is_vep_coding) + svannotate_hits
	if coding_hits:
		row_weights[("", "Coding")] = coding_hits

	if is_dbgap:
		row_weights[("", "dbGaP")] = 1
	if is_gnomad_matched:
		row_weights[("", "gnomAD Matched")] = 1

	return {row_key: weight for row_key, weight in row_weights.items() if weight > 0}


def get_genotype_weights(record):
	carrier_count = 0
	alt_allele_count = 0

	for sample in record.samples.values():
		genotype = sample.get("GT")
		if genotype is None:
			continue

		alt_alleles = sum(1 for allele in genotype if allele is not None and allele > 0)
		alt_allele_count += alt_alleles
		if alt_alleles > 0:
			carrier_count += 1

	return carrier_count, alt_allele_count


def write_table(path, table, integer_output):
	with open(path, "w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		writer.writerow(["annotation_group", "annotation"] + COLUMN_BUCKETS)
		for row_key in ROW_ORDER:
			values = []
			for column in COLUMN_BUCKETS:
				value = table[row_key][column]
				if integer_output:
					values.append(str(int(value)))
				else:
					values.append(str(value))
			writer.writerow([row_key[0], row_key[1]] + values)


site_table = init_table()
sample_table = init_table()
allele_table = init_table()

vcf_in = pysam.VariantFile(VCF_PATH)
consequence_idx = get_vep_consequence_index(vcf_in.header)
sample_count = len(vcf_in.header.samples)

with open(SAMPLE_COUNT_OUTPUT, "w") as handle:
	handle.write(f"{sample_count}\n")

for record in vcf_in:
	column = determine_column(record)
	row_weights = determine_row_weights(record, consequence_idx)

	if DO_PER_SAMPLE or DO_PER_ALLELE:
		carrier_count, alt_allele_count = get_genotype_weights(record)
	else:
		carrier_count, alt_allele_count = 0, 0

	for row_key, weight in row_weights.items():
		site_table[row_key][column] += weight
		if DO_PER_SAMPLE:
			sample_table[row_key][column] += carrier_count * weight
		if DO_PER_ALLELE:
			allele_table[row_key][column] += alt_allele_count * weight

write_table(SITE_OUTPUT, site_table, integer_output=True)
write_table(SAMPLE_OUTPUT, sample_table, integer_output=True)
write_table(ALLELE_OUTPUT, allele_table, integer_output=True)
PYCODE
	>>>

	output {
		File site_counts_tsv = "~{prefix}.sites.raw.tsv"
		File sample_counts_tsv = "~{prefix}.samples.raw.tsv"
		File allele_counts_tsv = "~{prefix}.alleles.raw.tsv"
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
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
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
OUTPUT = "~{prefix}.tsv"
PARENT_ROWS = {"ME", "DUP", "SVAnnotate"}


def format_label(key):
	annotation_group, annotation = key
	if annotation_group:
		return "", annotation
	if annotation in PARENT_ROWS:
		return annotation, "total"
	return annotation, ""


def format_site_value(value, total):
	count = int(round(value))
	percentage = 0.0 if total <= 0 else (value / total) * 100.0
	return str(count), f"({percentage:.2f}%)"


def format_normalized_value(value):
	return f"{value:.2f}"


header = None
row_order = []
counts = {}

for path in COUNT_FILES:
	with open(path, "r", newline="") as handle:
		reader = csv.reader(handle, delimiter="\t")
		current_header = next(reader)
		if header is None:
			header = current_header
		elif current_header != header:
			raise ValueError(f"Mismatched headers while merging count tables: {path}")

		for row in reader:
			key = (row[0], row[1])
			if key not in counts:
				counts[key] = [0.0] * (len(header) - 2)
				row_order.append(key)
			counts[key] = [left + float(right) for left, right in zip(counts[key], row[2:])]

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
	if MODE == "alleles":
		denominator *= 2.0

with open(OUTPUT, "w", newline="") as handle:
	writer = csv.writer(handle, delimiter="\t")
	writer.writerow(header)
	all_counts = counts.get(("", "All"))
	if all_counts is None:
		raise ValueError("Merged count tables are missing the All row required for output formatting")

	for key in row_order:
		values = counts[key]
		if MODE != "sites":
			values = [value / denominator for value in values]

		display_group, display_annotation = format_label(key)
		if MODE == "sites":
			formatted_pairs = [format_site_value(value, total) for value, total in zip(values, all_counts)]
			count_values = [count_value for count_value, _ in formatted_pairs]
			percentage_values = [percentage_value for _, percentage_value in formatted_pairs]
			writer.writerow([display_group, display_annotation] + count_values)
			writer.writerow(["", ""] + percentage_values)
		else:
			formatted_values = [format_normalized_value(value) for value in values]
			writer.writerow([display_group, display_annotation] + formatted_values)
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
