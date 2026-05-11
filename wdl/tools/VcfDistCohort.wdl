version 1.0

import "../utils/Helpers.wdl" as Helpers
import "../utils/Structs.wdl"
import "BackbonePhase.wdl" as BackbonePhase

workflow VcfDistCohort {
	input {
		Array[File] eval_vcfs
		Array[File] eval_vcf_idxs
		Array[File] truth_vcfs
		Array[File] truth_vcf_idxs
		Array[String] contigs
		String prefix

		File ref_fa

		Array[String]? subset_samples
		String? vcfdist_args

		String utils_docker
		String vcfdist_docker

		RuntimeAttr? runtime_attr_subset_assignment_samples
		RuntimeAttr? runtime_attr_assign_samples
		RuntimeAttr? runtime_attr_extract_assigned_samples
		RuntimeAttr? runtime_attr_subset_truth_contig
		RuntimeAttr? runtime_attr_subset_truth_samples
		RuntimeAttr? runtime_attr_subset_eval_samples
		RuntimeAttr? runtime_attr_vcfdist_shard
		RuntimeAttr? runtime_attr_aggregate_results
	}

	# Optionally subset eval VCF to a sample subset for assignment
	if (defined(subset_samples)) {
		call Helpers.SubsetVcfToSamples as SubsetAssignmentVcf {
			input:
				vcf = eval_vcfs[0],
				vcf_idx = eval_vcf_idxs[0],
				samples = select_first([subset_samples]),
				prefix = "~{prefix}.assignment_subset",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_subset_assignment_samples
		}
	}

	File assignment_vcf = select_first([SubsetAssignmentVcf.subset_vcf, eval_vcfs[0]])
	File assignment_vcf_idx = select_first([SubsetAssignmentVcf.subset_vcf_idx, eval_vcf_idxs[0]])

	# Assign each sample to its best-matching truth VCF
	call BackbonePhase.AssignSamplesToBaseVcfs {
		input:
			vcf = assignment_vcf,
			vcf_idx = assignment_vcf_idx,
			base_vcfs = truth_vcfs,
			prefix = "~{prefix}.assignments",
			docker = utils_docker,
			runtime_attr_override = runtime_attr_assign_samples
	}

	# Extract assigned sample lists per truth VCF
	scatter (j in range(length(truth_vcfs))) {
		call ExtractAssignedSamples {
			input:
				assignment_tsv = AssignSamplesToBaseVcfs.assignment_tsv,
				truth_vcf_index = j,
				prefix = "~{prefix}.truth_~{j}.assigned_samples",
				utils_docker = utils_docker,
				runtime_attr_override = runtime_attr_extract_assigned_samples
		}
	}

	# For each contig x truth VCF pair, subset to common samples and run vcfdist
	scatter (i in range(length(eval_vcfs))) {
		scatter (j in range(length(truth_vcfs))) {
			call Helpers.SubsetVcfToContig {
				input:
					vcf = truth_vcfs[j],
					vcf_idx = truth_vcf_idxs[j],
					contig = contigs[i],
					prefix = "~{prefix}.truth_~{j}.~{contigs[i]}.contig",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_truth_contig
			}

			call Helpers.SubsetVcfToSamples as SubsetTruthSamples {
				input:
					vcf = SubsetVcfToContig.subset_vcf,
					vcf_idx = SubsetVcfToContig.subset_vcf_idx,
					samples = ExtractAssignedSamples.assigned_samples[j],
					prefix = "~{prefix}.truth_~{j}.~{contigs[i]}.samples",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_truth_samples
			}

			call Helpers.SubsetVcfToSamples as SubsetEvalSamples {
				input:
					vcf = eval_vcfs[i],
					vcf_idx = eval_vcf_idxs[i],
					samples = ExtractAssignedSamples.assigned_samples[j],
					prefix = "~{prefix}.eval_~{contigs[i]}.truth_~{j}.samples",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_eval_samples
			}

			call RunVcfDistCohortShard {
				input:
					eval_vcf = SubsetEvalSamples.subset_vcf,
					eval_vcf_idx = SubsetEvalSamples.subset_vcf_idx,
					truth_vcf = SubsetTruthSamples.subset_vcf,
					truth_vcf_idx = SubsetTruthSamples.subset_vcf_idx,
					contig = contigs[i],
					ref_fa = ref_fa,
					vcfdist_args = select_first([vcfdist_args, ""]),
					prefix = "~{prefix}.~{contigs[i]}.truth_~{j}",
					docker = vcfdist_docker,
					runtime_attr_override = runtime_attr_vcfdist_shard
			}
		}

		Array[File] contig_phasing_tsvs = RunVcfDistCohortShard.phasing_tsv
		Array[File] contig_precision_recall_tsvs = RunVcfDistCohortShard.precision_recall_tsv
	}

	Array[File] all_phasing_tsvs = flatten(contig_phasing_tsvs)
	Array[File] all_precision_recall_tsvs = flatten(contig_precision_recall_tsvs)

	call AggregateVcfDistResults {
		input:
			phasing_tsvs = all_phasing_tsvs,
			precision_recall_tsvs = all_precision_recall_tsvs,
			contigs = contigs,
			prefix = prefix,
			utils_docker = utils_docker,
			runtime_attr_override = runtime_attr_aggregate_results
	}

	output {
		File phasing_table = AggregateVcfDistResults.phasing_table
		File precision_recall_table = AggregateVcfDistResults.precision_recall_table
		File phasing_plot = AggregateVcfDistResults.phasing_plot
		File missing_samples = AssignSamplesToBaseVcfs.missing_samples
	}
}

task ExtractAssignedSamples {
	input {
		File assignment_tsv
		Int truth_vcf_index
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
		if int(base_idx) == ~{truth_vcf_index}:
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

task RunVcfDistCohortShard {
	input {
		File eval_vcf
		File eval_vcf_idx
		File truth_vcf
		File truth_vcf_idx
		String contig
		File ref_fa
		String? vcfdist_args
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		ln -s "~{eval_vcf}" eval.vcf.gz
		ln -s "~{eval_vcf_idx}" eval.vcf.gz.tbi
		ln -s "~{truth_vcf}" truth.vcf.gz
		ln -s "~{truth_vcf_idx}" truth.vcf.gz.tbi
		ln -s "~{ref_fa}" ref.fa

		python3 <<CODE
import gzip
import subprocess
import os
import csv


def split_vcf(input_vcf, output_dir, suffix):
	"""Split a multi-sample VCF into per-sample single-sample VCFs in one pass."""
	os.makedirs(output_dir, exist_ok=True)
	meta_headers = []
	samples = []
	handles = {}
	sample_cols = {}

	with gzip.open(input_vcf, "rt") as f:
		for line in f:
			if line.startswith("##"):
				meta_headers.append(line)
			elif line.startswith("#CHROM"):
				fields = line.rstrip("\n").split("\t")
				samples = fields[9:]
				fixed = "\t".join(fields[:9])
				for i, sample in enumerate(samples):
					sample_cols[sample] = i + 9
					path = os.path.join(output_dir, f"{sample}.{suffix}.vcf")
					handle = open(path, "w")
					for mh in meta_headers:
						handle.write(mh)
					handle.write(fixed + "\t" + sample + "\n")
					handles[sample] = handle
			else:
				fields = line.rstrip("\n").split("\t")
				fixed = "\t".join(fields[:9])
				for sample, col in sample_cols.items():
					gt = fields[col].split(":")[0]
					alleles = gt.replace("|", "/").split("/")
					if all(a in ("0", ".") for a in alleles):
						continue
					handles[sample].write(fixed + "\t" + fields[col] + "\n")

	for h in handles.values():
		h.close()

	for sample in samples:
		path = os.path.join(output_dir, f"{sample}.{suffix}.vcf")
		subprocess.run(["bgzip", path], check=True)
		subprocess.run(["tabix", "-p", "vcf", f"{path}.gz"], check=True)

	return samples


# Split multi-sample VCFs into per-sample VCFs
eval_samples = split_vcf("eval.vcf.gz", "per_sample", "eval")
truth_samples = split_vcf("truth.vcf.gz", "per_sample", "truth")
common_samples = sorted(set(eval_samples) & set(truth_samples))

contig = "~{contig}"
vcfdist_extra = "~{if defined(vcfdist_args) then vcfdist_args else ""}"

phasing_rows = []
pr_rows = []

for sample in common_samples:
	eval_path = f"per_sample/{sample}.eval.vcf.gz"
	truth_path = f"per_sample/{sample}.truth.vcf.gz"
	result_dir = f"results/{sample}/"
	os.makedirs(result_dir, exist_ok=True)

	cmd = ["vcfdist", eval_path, truth_path, "ref.fa", "-p", result_dir]
	if vcfdist_extra:
		cmd.extend(vcfdist_extra.split())

	result = subprocess.run(cmd, capture_output=True, text=True)
	if result.returncode != 0:
		print(f"WARNING: vcfdist failed for {sample}: {result.stderr}", flush=True)
		continue

	# Collect phasing summary
	phasing_file = f"{result_dir}phasing-summary.tsv"
	if os.path.exists(phasing_file):
		with open(phasing_file) as f:
			reader = csv.DictReader(f, delimiter="\t")
			for row in reader:
				row["contig"] = contig
				row["sample"] = sample
				phasing_rows.append(dict(row))

	# Sum superclusters across phase blocks for rate computation
	phase_blocks_file = f"{result_dir}phase-blocks.tsv"
	total_sc = 0
	if os.path.exists(phase_blocks_file):
		with open(phase_blocks_file) as f:
			reader = csv.DictReader(f, delimiter="\t")
			for row in reader:
				total_sc += int(row.get("SUPERCLUSTERS", 0))
	if phasing_rows and phasing_rows[-1]["sample"] == sample:
		phasing_rows[-1]["TOTAL_SUPERCLUSTERS"] = str(total_sc)

	# Collect precision-recall summary
	pr_file = f"{result_dir}precision-recall-summary.tsv"
	if os.path.exists(pr_file):
		with open(pr_file) as f:
			reader = csv.DictReader(f, delimiter="\t")
			for row in reader:
				row["contig"] = contig
				row["sample"] = sample
				pr_rows.append(dict(row))

prefix = "~{prefix}"

# Write shard phasing TSV
with open(f"{prefix}.phasing.tsv", "w") as f:
	if phasing_rows:
		base_cols = [k for k in phasing_rows[0] if k not in ("contig", "sample")]
		cols = ["contig", "sample"] + base_cols
		f.write("\t".join(cols) + "\n")
		for row in phasing_rows:
			f.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
	else:
		f.write("contig\tsample\n")

# Write shard precision-recall TSV
with open(f"{prefix}.precision_recall.tsv", "w") as f:
	if pr_rows:
		base_cols = [k for k in pr_rows[0] if k not in ("contig", "sample")]
		cols = ["contig", "sample"] + base_cols
		f.write("\t".join(cols) + "\n")
		for row in pr_rows:
			f.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
	else:
		f.write("contig\tsample\n")
CODE
	>>>

	output {
		File phasing_tsv = "~{prefix}.phasing.tsv"
		File precision_recall_tsv = "~{prefix}.precision_recall.tsv"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 4,
		mem_gb: 32,
		disk_gb: 4 * (ceil(size(eval_vcf, "GiB")) + ceil(size(truth_vcf, "GiB"))) + ceil(size(ref_fa, "GiB")) + 50,
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

task AggregateVcfDistResults {
	input {
		Array[File] phasing_tsvs
		Array[File] precision_recall_tsvs
		Array[String] contigs
		String prefix
		String utils_docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<CODE
from collections import defaultdict
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


with open("~{write_lines(contigs)}") as handle:
	ordered_contigs = sorted({line.strip() for line in handle if line.strip()}, key=contig_sort_key)

# Load and deduplicate phasing rows across shards
phasing_rows = []
phasing_header = None
with open("~{write_lines(phasing_tsvs)}") as manifest:
	for path in manifest:
		path = path.strip()
		if not path:
			continue
		with open(path) as f:
			reader = csv.DictReader(f, delimiter="\t")
			if phasing_header is None and reader.fieldnames and len(reader.fieldnames) > 2:
				phasing_header = reader.fieldnames
			for row in reader:
				phasing_rows.append(row)

seen_phasing = set()
unique_phasing = []
for row in phasing_rows:
	key = (row.get("contig", ""), row.get("sample", ""))
	if key not in seen_phasing:
		seen_phasing.add(key)
		unique_phasing.append(row)
unique_phasing.sort(key=lambda r: (contig_sort_key(r.get("contig", "")), r.get("sample", "")))

with open("~{prefix}.phasing.tsv", "w") as f:
	if phasing_header and unique_phasing:
		f.write("\t".join(phasing_header) + "\n")
		for row in unique_phasing:
			f.write("\t".join(str(row.get(c, "")) for c in phasing_header) + "\n")
	else:
		f.write("contig\tsample\n")

# Load and deduplicate precision-recall rows
pr_rows = []
pr_header = None
with open("~{write_lines(precision_recall_tsvs)}") as manifest:
	for path in manifest:
		path = path.strip()
		if not path:
			continue
		with open(path) as f:
			reader = csv.DictReader(f, delimiter="\t")
			if pr_header is None and reader.fieldnames and len(reader.fieldnames) > 2:
				pr_header = reader.fieldnames
			for row in reader:
				pr_rows.append(row)

seen_pr = set()
unique_pr = []
for row in pr_rows:
	key = (row.get("contig", ""), row.get("sample", ""), row.get("VAR_TYPE", ""), row.get("MIN_QUAL", ""))
	if key not in seen_pr:
		seen_pr.add(key)
		unique_pr.append(row)
unique_pr.sort(key=lambda r: (contig_sort_key(r.get("contig", "")), r.get("sample", ""), r.get("VAR_TYPE", "")))

with open("~{prefix}.precision_recall.tsv", "w") as f:
	if pr_header and unique_pr:
		f.write("\t".join(pr_header) + "\n")
		for row in unique_pr:
			f.write("\t".join(str(row.get(c, "")) for c in pr_header) + "\n")
	else:
		f.write("contig\tsample\n")

# SVG plot of switch error rate per chromosome
grouped_rates = {contig: [] for contig in ordered_contigs}
for row in unique_phasing:
	contig = row.get("contig", "")
	switch_errors = int(row.get("SWITCH_ERRORS", 0))
	total_sc = int(row.get("TOTAL_SUPERCLUSTERS", 0))
	if total_sc > 0:
		rate_pct = 100.0 * switch_errors / float(total_sc)
		grouped_rates[contig].append(rate_pct)

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


title = "Switch Error Rate (vcfdist)"
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
	parts.append(f'<text x="{width / 2:.2f}" y="{height / 2:.2f}" font-size="14" text-anchor="middle" font-family="Arial">No phased superclusters found for comparison</text>')

parts.append("</svg>")
with open("~{prefix}.switch_error_rate.svg", "w") as out:
	out.write("\n".join(parts))
CODE
	>>>

	output {
		File phasing_table = "~{prefix}.phasing.tsv"
		File precision_recall_table = "~{prefix}.precision_recall.tsv"
		File phasing_plot = "~{prefix}.switch_error_rate.svg"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: ceil(size(phasing_tsvs, "GiB")) + ceil(size(precision_recall_tsvs, "GiB")) + 15,
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
