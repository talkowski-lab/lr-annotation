version 1.0

import "../utils/Helpers.wdl" as Helpers
import "../tools/BackbonePhase.wdl" as BackbonePhase

workflow PlotPhasingResults {
	input {
		Array[File] backbone_phased_vcfs
		Array[File] backbone_phased_vcf_idxs
		Array[File] base_vcfs
		Array[File] base_vcf_idxs
		Array[String] contigs
		String prefix

		Int max_variants = -1
		Array[String]? subset_samples

		String utils_docker

		RuntimeAttr? runtime_attr_subset_assignment_samples
		RuntimeAttr? runtime_attr_assign_samples
		RuntimeAttr? runtime_attr_extract_assigned_samples
		RuntimeAttr? runtime_attr_subset_base_contig
		RuntimeAttr? runtime_attr_subset_base_samples
		RuntimeAttr? runtime_attr_subset_backbone_samples
		RuntimeAttr? runtime_attr_compare_shard
		RuntimeAttr? runtime_attr_build_vcf_table
		RuntimeAttr? runtime_attr_aggregate_results
	}

	if (defined(subset_samples)) {
		call Helpers.SubsetVcfToSamples as SubsetAssignmentVcf {
			input:
				vcf = backbone_phased_vcfs[0],
				vcf_idx = backbone_phased_vcf_idxs[0],
				samples = select_first([subset_samples]),
				prefix = "~{prefix}.assignment_subset",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_subset_assignment_samples
		}
	}

	File assignment_vcf = select_first([SubsetAssignmentVcf.subset_vcf, backbone_phased_vcfs[0]])
	File assignment_vcf_idx = select_first([SubsetAssignmentVcf.subset_vcf_idx, backbone_phased_vcf_idxs[0]])

	call BackbonePhase.AssignSamplesToBaseVcfs as AssignSamplesToBaseVcfs {
		input:
			vcf = assignment_vcf,
			vcf_idx = assignment_vcf_idx,
			base_vcfs = base_vcfs,
			prefix = "~{prefix}.assignments",
			docker = utils_docker,
			runtime_attr_override = runtime_attr_assign_samples
	}

	scatter (j in range(length(base_vcfs))) {
		call ExtractAssignedSamplesForBaseVcf {
			input:
				assignment_tsv = AssignSamplesToBaseVcfs.assignment_tsv,
				base_vcf_index = j,
				prefix = "~{prefix}.base_~{j}.assigned_samples",
				utils_docker = utils_docker,
				runtime_attr_override = runtime_attr_extract_assigned_samples
		}
	}

	scatter (i in range(length(backbone_phased_vcfs))) {
		scatter (j in range(length(base_vcfs))) {
			call Helpers.SubsetVcfToContig {
				input:
					vcf = base_vcfs[j],
					vcf_idx = base_vcf_idxs[j],
					contig = contigs[i],
					prefix = "~{prefix}.base_~{j}.~{contigs[i]}.contig",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_base_contig
			}

			call Helpers.SubsetVcfToSamples as SubsetBaseSamples {
				input:
					vcf = SubsetVcfToContig.subset_vcf,
					vcf_idx = SubsetVcfToContig.subset_vcf_idx,
					samples = ExtractAssignedSamplesForBaseVcf.assigned_samples[j],
					prefix = "~{prefix}.base_~{j}.~{contigs[i]}.samples",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_base_samples
			}

			call Helpers.SubsetVcfToSamples as SubsetBackboneSamples {
				input:
					vcf = backbone_phased_vcfs[i],
					vcf_idx = backbone_phased_vcf_idxs[i],
					samples = ExtractAssignedSamplesForBaseVcf.assigned_samples[j],
					prefix = "~{prefix}.backbone_~{contigs[i]}.base_~{j}.samples",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_backbone_samples
			}

			call CompareBackbonePhasingShard {
				input:
					backbone_phased_vcf = SubsetBackboneSamples.subset_vcf,
					backbone_phased_vcf_idx = SubsetBackboneSamples.subset_vcf_idx,
					base_vcf = SubsetBaseSamples.subset_vcf,
					base_vcf_idx = SubsetBaseSamples.subset_vcf_idx,
					contig = contigs[i],
					max_variants = max_variants,
					prefix = "~{prefix}.~{contigs[i]}.base_~{j}",
					utils_docker = utils_docker,
					runtime_attr_override = runtime_attr_compare_shard
			}
		}

		Array[File] contig_outside_tr_tsvs = CompareBackbonePhasingShard.outside_tr_tsv
		Array[File] contig_tr_enveloped_tsvs = CompareBackbonePhasingShard.tr_enveloped_tsv

		call BuildContigVcfTable {
			input:
				backbone_phased_vcf = backbone_phased_vcfs[i],
				backbone_phased_vcf_idx = backbone_phased_vcf_idxs[i],
				status_tsv_gzs = CompareBackbonePhasingShard.status_tsv_gz,
				missing_samples_file = AssignSamplesToBaseVcfs.missing_samples,
				subset_samples = subset_samples,
				prefix = "~{prefix}.~{contigs[i]}.variant_status",
				utils_docker = utils_docker,
				runtime_attr_override = runtime_attr_build_vcf_table
		}
	}

	Array[File] outside_tr_tsvs = flatten(contig_outside_tr_tsvs)
	Array[File] tr_enveloped_tsvs = flatten(contig_tr_enveloped_tsvs)
	Array[File] variant_status_tables = BuildContigVcfTable.vcf_table_tsv_gz

	call AggregatePhasingResults {
		input:
			outside_tr_tsvs = outside_tr_tsvs,
			tr_enveloped_tsvs = tr_enveloped_tsvs,
			contigs = contigs,
			prefix = prefix,
			utils_docker = utils_docker,
			runtime_attr_override = runtime_attr_aggregate_results
	}

	output {
		File outside_tr_table = AggregatePhasingResults.outside_tr_table
		File tr_enveloped_table = AggregatePhasingResults.tr_enveloped_table
		File outside_tr_plot = AggregatePhasingResults.outside_tr_plot
		File tr_enveloped_plot = AggregatePhasingResults.tr_enveloped_plot
		File missing_samples = AssignSamplesToBaseVcfs.missing_samples
		Array[File] vcf_tables = variant_status_tables
	}
}

task ExtractAssignedSamplesForBaseVcf {
	input {
		File assignment_tsv
		Int base_vcf_index
		String prefix
		String utils_docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<CODE
samples = []
with open("~{assignment_tsv}") as handle:
	for line in handle:
		sample, base_idx = line.rstrip("\n").split("\t")
		if int(base_idx) == ~{base_vcf_index}:
			samples.append(sample)

samples = sorted(samples)
with open("~{prefix}.samples.txt", "w") as out:
	for sample in samples:
		out.write(sample + "\n")
CODE
	>>>

	output {
		File samples_file = "~{prefix}.samples.txt"
		Array[String] assigned_samples = read_lines("~{prefix}.samples.txt")
		Int sample_count = length(read_lines("~{prefix}.samples.txt"))
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: 10,
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
		docker: utils_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task CompareBackbonePhasingShard {
	input {
		File backbone_phased_vcf
		File backbone_phased_vcf_idx
		File base_vcf
		File base_vcf_idx
		String contig
		Int max_variants
		String prefix
		String utils_docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		ln -s "~{backbone_phased_vcf}" backbone.vcf.gz
		ln -s "~{backbone_phased_vcf_idx}" backbone.vcf.gz.tbi
		ln -s "~{base_vcf}" base.vcf.gz
		ln -s "~{base_vcf_idx}" base.vcf.gz.tbi

		python3 <<CODE
from collections import defaultdict
import gzip
import pysam


def get_allele_type(record):
	try:
		return record.info["allele_type"]
	except (KeyError, ValueError):
		return ""


def parse_sample_gt(sample_data, require_phased=False):
	gt = sample_data.get("GT")
	if gt is None or None in gt or len(gt) != 2:
		return None
	if gt[0] == gt[1]:
		return None
	if require_phased and not sample_data.phased:
		return None
	return tuple(int(allele) for allele in gt), sample_data.phased


def normalize_biallelic_variant(pos, ref, alt):
	pos = int(pos)
	ref = ref.upper()
	alt = alt.upper()
	while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
		ref = ref[1:]
		alt = alt[1:]
		pos += 1
	while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
		ref = ref[:-1]
		alt = alt[:-1]
	return pos, ref, alt


def iter_comparable_calls(record, sample_name, require_phased=False):
	parsed = parse_sample_gt(record.samples[sample_name], require_phased=require_phased)
	if parsed is None:
		return []
	gt, _phased = parsed
	allele_type = get_allele_type(record)
	alt_values = [alt.upper() for alt in record.alts] if record.alts else []
	if allele_type == "trv":
		calls = []
		for alt_index, alt in enumerate(alt_values, start=1):
			biallelic_gt = tuple(1 if allele == alt_index else 0 for allele in gt)
			if biallelic_gt[0] == biallelic_gt[1]:
				continue
			norm_pos, norm_ref, norm_alt = normalize_biallelic_variant(record.pos, record.ref, alt)
			calls.append({
				"pos": norm_pos,
				"gt": biallelic_gt,
				"key": (record.contig, norm_pos, norm_ref, (norm_alt,)),
			})
		return calls

	return [{
		"pos": record.pos,
		"gt": gt,
		"key": (record.contig, record.pos, record.ref.upper(), tuple(alt_values)),
	}]


def normalize_record_key(record):
	alt_values = [alt.upper() for alt in record.alts] if record.alts else []
	allele_type = get_allele_type(record)
	return (
		record.contig,
		str(record.pos),
		record.id if record.id else ".",
		record.ref.upper(),
		",".join(alt_values),
		"" if allele_type in ("", ".") else allele_type,
	)


def summarize(points):
	points.sort(key=lambda item: item[0])
	states = [state for _pos, state, _record_key in points]
	matched_count = len(states)
	if matched_count == 0:
		return 0, 0, {}
	if sum(states) < matched_count / 2.0:
		states = [1 - state for state in states]
	status_priority = {"XC": 0, "CN": 1, "SW": 2}
	record_statuses = {}
	switch_error_count = 0
	for idx, (_pos, _state, record_key) in enumerate(points):
		if idx > 0 and states[idx] != states[idx - 1]:
			status = "SW"
			switch_error_count += 1
		elif states[idx] == 1:
			status = "XC"
		else:
			status = "CN"
		current = record_statuses.get(record_key)
		if current is None or status_priority[status] > status_priority[current]:
			record_statuses[record_key] = status
	return matched_count, switch_error_count, record_statuses


with pysam.VariantFile("backbone.vcf.gz") as backbone_in:
	samples = list(backbone_in.header.samples)

base_calls = defaultdict(dict)
with pysam.VariantFile("base.vcf.gz") as base_in:
	base_sample_set = set(base_in.header.samples)
	for record in base_in:
		if not record.alts:
			continue
		for sample in samples:
			if sample not in base_sample_set:
				continue
			for call in iter_comparable_calls(record, sample, require_phased=True):
				base_calls[sample][call["key"]] = call["gt"]

collection_points = {
	"outside_tr": defaultdict(list),
	"tr_enveloped": defaultdict(list),
}

limit = ~{max_variants}
seen_records = 0
with pysam.VariantFile("backbone.vcf.gz") as backbone_in:
	for record in backbone_in:
		if not record.alts:
			continue
		if limit >= 0 and seen_records >= limit:
			break
		seen_records += 1
		allele_type = get_allele_type(record)
		collection = "tr_enveloped" if allele_type == "trv" else "outside_tr"
		normalized_record_key = normalize_record_key(record)
		for sample in samples:
			for call in iter_comparable_calls(record, sample, require_phased=True):
				base_gt = base_calls[sample].get(call["key"])
				if base_gt is None:
					continue
				if call["gt"] == base_gt:
					state = 1
				elif call["gt"] == (base_gt[1], base_gt[0]):
					state = 0
				else:
					continue
				collection_points[collection][sample].append((call["pos"], state, normalized_record_key))

with open("~{prefix}.outside_tr.tsv", "w") as outside_out, open("~{prefix}.tr_enveloped.tsv", "w") as tr_out, gzip.open("~{prefix}.status.tsv.gz", "wt") as status_out:
	outside_out.write("contig\tsample\tmatched_count\tswitch_error_count\n")
	tr_out.write("contig\tsample\tmatched_count\tswitch_error_count\n")
	status_out.write("chrom\tpos\tvariant_id\tref\talt\tallele_type\tsample\tstatus\n")
	for collection, handle in (("outside_tr", outside_out), ("tr_enveloped", tr_out)):
		for sample in samples:
			matched_count, switch_error_count, record_statuses = summarize(collection_points[collection][sample])
			handle.write(f"~{contig}\t{sample}\t{matched_count}\t{switch_error_count}\n")
			for record_key, status in record_statuses.items():
				status_out.write("\t".join((*record_key, sample, status)) + "\n")
CODE
	>>>

	output {
		File outside_tr_tsv = "~{prefix}.outside_tr.tsv"
		File tr_enveloped_tsv = "~{prefix}.tr_enveloped.tsv"
		File status_tsv_gz = "~{prefix}.status.tsv.gz"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 12,
		disk_gb: 50 * ceil(size(backbone_phased_vcf, "GiB")) + ceil(size(base_vcf, "GiB")) + 25,
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
		docker: utils_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task BuildContigVcfTable {
	input {
		File backbone_phased_vcf
		File backbone_phased_vcf_idx
		Array[File] status_tsv_gzs
		File missing_samples_file
		Array[String]? subset_samples
		String prefix
		String utils_docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		ln -s "~{backbone_phased_vcf}" backbone.vcf.gz
		ln -s "~{backbone_phased_vcf_idx}" backbone.vcf.gz.tbi

		python3 <<CODE
import gzip
import pysam


def get_allele_type(rec):
	try:
		return rec.info["allele_type"]
	except (KeyError, ValueError):
		return ""


def make_record_key(rec):
	alt_values = [alt.upper() for alt in rec.alts] if rec.alts else []
	return (
		rec.contig,
		str(rec.pos),
		rec.id if rec.id else ".",
		rec.ref.upper(),
		",".join(alt_values),
		get_allele_type(rec),
	)


def classify_backbone_status(sample_data, comparable_status, is_missing_sample):
	gt = sample_data.get("GT")
	if gt is None:
		return "NC"
	if len(gt) != 2 or any(allele is None for allele in gt):
		return "NC"
	if gt[0] == gt[1]:
		if gt[0] == 0:
			return "."
		if gt[0] > 0:
			return "HA"
		return "OTH"
	if not sample_data.phased:
		return "UP"
	if comparable_status is not None:
		return comparable_status
	if is_missing_sample:
		return "MS"
	if gt[0] != gt[1]:
		return "NF"
	return "OTH"


subset_samples = {line for line in """~{sep='\n' select_first([subset_samples, []])}""".splitlines() if line}
missing_samples = {line.strip() for line in open("~{missing_samples_file}") if line.strip()}

status_priority = {"XC": 0, "CN": 1, "SW": 2}
status_by_key = {}
with open("~{write_lines(status_tsv_gzs)}") as manifest:
	for path in manifest:
		path = path.strip()
		if not path:
			continue
		with gzip.open(path, "rt") as handle:
			next(handle)
			for line in handle:
				chrom, pos, variant_id, ref, alt, allele_type, sample, status = line.rstrip("\n").split("\t")
				key = ((chrom, pos, variant_id, ref, alt, allele_type), sample)
				current = status_by_key.get(key)
				if current is None or status_priority[status] > status_priority[current]:
					status_by_key[key] = status

with pysam.VariantFile("backbone.vcf.gz") as backbone_in:
	all_samples = list(backbone_in.header.samples)
	samples = [sample for sample in all_samples if not subset_samples or sample in subset_samples]
	with gzip.open("~{prefix}.tsv.gz", "wt") as out:
		out.write("chrom\tpos\tvariant_id\tref\talt\tallele_type\t" + "\t".join(samples) + "\n")
		for rec in backbone_in:
			record_key = make_record_key(rec)
			allele_type = record_key[5] if record_key[5] else "."
			row = [record_key[0], record_key[1], record_key[2], record_key[3], record_key[4], allele_type]
			for sample in samples:
				status = classify_backbone_status(
					rec.samples[sample],
					status_by_key.get((record_key, sample)),
					sample in missing_samples,
				)
				row.append(status)
			out.write("\t".join(row) + "\n")
CODE
	>>>

	output {
		File vcf_table_tsv_gz = "~{prefix}.tsv.gz"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 6,
		disk_gb: ceil(size(backbone_phased_vcf, "GiB")) + ceil(size(status_tsv_gzs, "GiB")) + 20,
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
		docker: utils_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task AggregatePhasingResults {
	input {
		Array[File] outside_tr_tsvs
		Array[File] tr_enveloped_tsvs
		Array[String] contigs
		String prefix
		String utils_docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<CODE
from html import escape
import csv
import math
import random


def contig_sort_key(contig):
	value = contig[3:] if contig.lower().startswith("chr") else contig
	lower_value = value.lower()
	if lower_value.isdigit():
		return (0, int(lower_value), contig)
	if lower_value == "x":
		return (1, 23, contig)
	if lower_value == "y":
		return (1, 24, contig)
	if lower_value in ("m", "mt"):
		return (1, 25, contig)
	return (2, lower_value, contig)


def load_rows(paths):
	rows_by_key = {}
	for path in paths:
		with open(path) as handle:
			reader = csv.DictReader(handle, delimiter="\t")
			for row in reader:
				key = (row["contig"], row["sample"])
				rows_by_key[key] = {
					"contig": row["contig"],
					"sample": row["sample"],
					"matched_count": rows_by_key.get(key, {}).get("matched_count", 0) + int(row["matched_count"]),
					"switch_error_count": rows_by_key.get(key, {}).get("switch_error_count", 0) + int(row["switch_error_count"]),
				}
	return list(rows_by_key.values())


def write_rows(rows, out_path):
	rows.sort(key=lambda row: (contig_sort_key(row["contig"]), row["sample"]))
	with open(out_path, "w") as out:
		out.write("contig\tsample\tmatched_count\tswitch_error_count\n")
		for row in rows:
			out.write(
				f"{row['contig']}\t{row['sample']}\t{row['matched_count']}\t{row['switch_error_count']}\n"
			)


def write_svg(rows, ordered_contigs, title, out_path):
	rng = random.Random(0)
	width = max(900, 140 * max(1, len(ordered_contigs)))
	height = 520
	margin_left = 80
	margin_right = 30
	margin_top = 55
	margin_bottom = 125
	plot_width = width - margin_left - margin_right
	plot_height = height - margin_top - margin_bottom
	plot_bottom = height - margin_bottom
	plot_top = margin_top

	grouped_rates = {contig: [] for contig in ordered_contigs}
	for row in rows:
		if row["matched_count"] <= 0:
			continue
		rate_pct = 100.0 * row["switch_error_count"] / float(row["matched_count"])
		grouped_rates[row["contig"]].append(rate_pct)

	all_rates = [rate for rates in grouped_rates.values() for rate in rates]
	positive_rates = [rate for rate in all_rates if rate > 0]

	def build_log_ticks(min_rate, max_rate):
		if max_rate <= 0:
			return [0.1, 0.2, 0.5, 1.0]
		lower_exp = math.floor(math.log10(min_rate))
		upper_exp = math.ceil(math.log10(max_rate))
		candidates = []
		for exp in range(lower_exp - 1, upper_exp + 2):
			base = 10 ** exp
			for multiplier in (1, 2, 5):
				value = multiplier * base
				if min_rate <= value <= max_rate:
					candidates.append(value)
		if not candidates:
			candidates = [min_rate, max_rate]
		return sorted(set(candidates))

	if positive_rates:
		min_positive_rate = min(positive_rates)
		max_rate = max(positive_rates)
		y_min = 10 ** math.floor(math.log10(min_positive_rate))
		y_max = 10 ** math.ceil(math.log10(max_rate))
		if y_min == y_max:
			y_min /= 10.0
			y_max *= 10.0
		ticks = build_log_ticks(y_min, y_max)
	else:
		y_min = 0.1
		y_max = 1.0
		ticks = [0.1, 0.2, 0.5, 1.0]

	def x_center(idx):
		if len(ordered_contigs) == 1:
			return margin_left + plot_width / 2.0
		return margin_left + idx * plot_width / float(len(ordered_contigs) - 1)

	def y_coord(rate_pct):
		clamped_rate = min(max(rate_pct, y_min), y_max)
		log_span = math.log10(y_max) - math.log10(y_min)
		return plot_bottom - ((math.log10(clamped_rate) - math.log10(y_min)) / log_span) * plot_height

	parts = [
		f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
		'<rect width="100%" height="100%" fill="white"/>',
		f'<text x="{width / 2:.2f}" y="28" font-size="18" text-anchor="middle" font-family="Arial">{escape(title)}</text>',
	]

	for tick in ticks:
		y = y_coord(tick)
		parts.append(f'<line x1="{margin_left}" y1="{y:.2f}" x2="{width - margin_right}" y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>')
		parts.append(f'<text x="{margin_left - 10}" y="{y + 4:.2f}" font-size="11" text-anchor="end" font-family="Arial">{tick:g}%</text>')

	parts.append(f'<line x1="{margin_left}" y1="{plot_top}" x2="{margin_left}" y2="{plot_bottom}" stroke="#222222" stroke-width="1.5"/>')
	parts.append(f'<line x1="{margin_left}" y1="{plot_bottom}" x2="{width - margin_right}" y2="{plot_bottom}" stroke="#222222" stroke-width="1.5"/>')
	parts.append(f'<text x="{width / 2:.2f}" y="{height - 25}" font-size="14" text-anchor="middle" font-family="Arial">Chromosome</text>')
	parts.append(f'<text x="24" y="{height / 2:.2f}" font-size="14" text-anchor="middle" transform="rotate(-90 24 {height / 2:.2f})" font-family="Arial">Switch Error Rate</text>')

	for idx, contig in enumerate(ordered_contigs):
		x = x_center(idx)
		parts.append(f'<text x="{x:.2f}" y="{plot_bottom + 18}" font-size="11" text-anchor="middle" font-family="Arial">{escape(contig)}</text>')
		rates = grouped_rates[contig]
		for rate_pct in rates:
			jitter = rng.uniform(-0.18, 0.18)
			jitter_span = plot_width / float(max(2, len(ordered_contigs)))
			cx = x + jitter * jitter_span
			cy = y_coord(rate_pct)
			parts.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="4" fill="#4c78a8" fill-opacity="0.55"/>')
		if rates:
			mean_rate = sum(rates) / float(len(rates))
			cy = y_coord(mean_rate)
			diamond = [
				(f"{x:.2f}", f"{cy - 7:.2f}"),
				(f"{x + 7:.2f}", f"{cy:.2f}"),
				(f"{x:.2f}", f"{cy + 7:.2f}"),
				(f"{x - 7:.2f}", f"{cy:.2f}"),
			]
			parts.append('<polygon points="' + ' '.join(','.join(point) for point in diamond) + '" fill="#c44e52"/>')
			parts.append(f'<text x="{x:.2f}" y="{plot_bottom + 34:.2f}" font-size="11" text-anchor="middle" fill="#c44e52" font-family="Arial">{mean_rate:.2f}%</text>')

	if not all_rates:
		parts.append(f'<text x="{width / 2:.2f}" y="{height / 2:.2f}" font-size="14" text-anchor="middle" font-family="Arial">No comparable phased het genotypes found</text>')

	parts.append('</svg>')
	with open(out_path, "w") as out:
		out.write("\n".join(parts))


with open("~{write_lines(contigs)}") as handle:
	ordered_contigs = sorted({line.strip() for line in handle if line.strip()}, key=contig_sort_key)

with open("~{write_lines(outside_tr_tsvs)}") as handle:
	outside_paths = [line.strip() for line in handle if line.strip()]

with open("~{write_lines(tr_enveloped_tsvs)}") as handle:
	tr_enveloped_paths = [line.strip() for line in handle if line.strip()]

outside_rows = load_rows(outside_paths)
tr_enveloped_rows = load_rows(tr_enveloped_paths)

write_rows(outside_rows, "~{prefix}.outside_tr.tsv")
write_rows(tr_enveloped_rows, "~{prefix}.tr_enveloped.tsv")
write_svg(outside_rows, ordered_contigs, "Outside TR", "~{prefix}.outside_tr.svg")
write_svg(tr_enveloped_rows, ordered_contigs, "TR Enveloped", "~{prefix}.tr_enveloped.svg")
CODE
	>>>

	output {
		File outside_tr_table = "~{prefix}.outside_tr.tsv"
		File tr_enveloped_table = "~{prefix}.tr_enveloped.tsv"
		File outside_tr_plot = "~{prefix}.outside_tr.svg"
		File tr_enveloped_plot = "~{prefix}.tr_enveloped.svg"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: ceil(size(outside_tr_tsvs, "GiB")) + ceil(size(tr_enveloped_tsvs, "GiB")) + 15,
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
		docker: utils_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}
