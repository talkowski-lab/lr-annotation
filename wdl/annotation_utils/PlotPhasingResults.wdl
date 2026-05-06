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
		String utils_docker

		Array[String]? subset_samples
		Int? max_variants

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
		Int? max_variants
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

		bcftools query -l backbone.vcf.gz > backbone.samples.txt
		bcftools query -f '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/allele_type\t%SAMPLE\t%GT\n]' backbone.vcf.gz > backbone.query.tsv
		bcftools query -f '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\n]' base.vcf.gz > base.query.tsv

		python3 <<CODE
from collections import defaultdict
import gzip


def parse_sample_gt(gt_string, require_phased=False):
	sep = "|" if "|" in gt_string else "/"
	parts = gt_string.split(sep)
	if len(parts) != 2:
		return None
	if "." in parts:
		return None
	try:
		alleles = tuple(int(part) for part in parts)
	except ValueError:
		return None
	if alleles[0] == alleles[1]:
		return None
	phased = sep == "|"
	if require_phased and not phased:
		return None
	return alleles, phased


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


def iter_comparable_calls(chrom, pos, ref, alt_string, allele_type, gt_string, require_phased=False):
	parsed = parse_sample_gt(gt_string, require_phased=require_phased)
	if parsed is None:
		return []
	gt, _phased = parsed
	allele_type = "" if allele_type in ("", ".") else allele_type
	alt_values = [alt.upper() for alt in alt_string.split(",") if alt]
	if allele_type == "trv":
		calls = []
		for alt_index, alt in enumerate(alt_values, start=1):
			biallelic_gt = tuple(1 if allele == alt_index else 0 for allele in gt)
			if biallelic_gt[0] == biallelic_gt[1]:
				continue
			norm_pos, norm_ref, norm_alt = normalize_biallelic_variant(pos, ref, alt)
			calls.append({
				"pos": norm_pos,
				"gt": biallelic_gt,
				"key": (chrom, norm_pos, norm_ref, (norm_alt,)),
			})
		return calls

	return [{
		"pos": int(pos),
		"gt": gt,
		"key": (chrom, int(pos), ref.upper(), tuple(alt_values)),
	}]


def normalize_record_key(chrom, pos, variant_id, ref, alt_string, allele_type):
	alt_values = [alt.upper() for alt in alt_string.split(",") if alt]
	return (
		chrom,
		str(pos),
		variant_id if variant_id else ".",
		ref.upper(),
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


samples = [line.strip() for line in open("backbone.samples.txt") if line.strip()]
base_calls = defaultdict(dict)
for line in open("base.query.tsv"):
	chrom, pos, _variant_id, ref, alt_string, sample, gt_string = line.rstrip("\n").split("\t")
	for call in iter_comparable_calls(chrom, pos, ref, alt_string, "", gt_string, require_phased=True):
		base_calls[sample][call["key"]] = call["gt"]

collection_points = {
	"outside_tr": defaultdict(list),
	"tr_enveloped": defaultdict(list),
}

limit = ~{if defined(max_variants) then max_variants else -1}
seen_records = 0
last_record = None

for line in open("backbone.query.tsv"):
	chrom, pos, variant_id, ref, alt_string, allele_type, sample, gt_string = line.rstrip("\n").split("\t")
	record_key = (chrom, pos, variant_id, ref, alt_string, allele_type)
	if record_key != last_record:
		if limit >= 0 and seen_records >= limit:
			break
		seen_records += 1
		last_record = record_key
	collection = "tr_enveloped" if allele_type == "trv" else "outside_tr"
	normalized_record_key = normalize_record_key(chrom, pos, variant_id, ref, alt_string, allele_type)
	for call in iter_comparable_calls(chrom, pos, ref, alt_string, allele_type, gt_string, require_phased=True):
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

status_rows = []
for collection in ("outside_tr", "tr_enveloped"):
	with open(f"~{prefix}.{collection}.tsv", "w") as out:
		out.write("contig\tsample\tmatched_count\tswitch_error_count\n")
		for sample in samples:
			matched_count, switch_error_count, record_statuses = summarize(collection_points[collection][sample])
			out.write(f"~{contig}\t{sample}\t{matched_count}\t{switch_error_count}\n")
			for record_key, status in record_statuses.items():
				status_rows.append((*record_key, sample, status))

with gzip.open("~{prefix}.status.tsv.gz", "wt") as out:
	out.write("chrom\tpos\tvariant_id\tref\talt\tallele_type\tsample\tstatus\n")
	for row in status_rows:
		out.write("\t".join(row) + "\n")
CODE
	>>>

	output {
		File outside_tr_tsv = "~{prefix}.outside_tr.tsv"
		File tr_enveloped_tsv = "~{prefix}.tr_enveloped.tsv"
		File status_tsv_gz = "~{prefix}.status.tsv.gz"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 5 * ceil(size(backbone_phased_vcf, "GiB")) + ceil(size(base_vcf, "GiB")) + 20,
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
	max_rate = max(all_rates) if all_rates else 0.0
	y_max = max(1.0, max_rate * 1.15)
	tick_count = 5
	ticks = [y_max * idx / float(tick_count) for idx in range(tick_count + 1)]

	def x_center(idx):
		if len(ordered_contigs) == 1:
			return margin_left + plot_width / 2.0
		return margin_left + idx * plot_width / float(len(ordered_contigs) - 1)

	def y_coord(rate_pct):
		return plot_bottom - (rate_pct / y_max) * plot_height

	parts = [
		f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
		'<rect width="100%" height="100%" fill="white"/>',
		f'<text x="{width / 2:.2f}" y="28" font-size="18" text-anchor="middle" font-family="Arial">{escape(title)}</text>',
	]

	for tick in ticks:
		y = y_coord(tick)
		parts.append(f'<line x1="{margin_left}" y1="{y:.2f}" x2="{width - margin_right}" y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>')
		parts.append(f'<text x="{margin_left - 10}" y="{y + 4:.2f}" font-size="11" text-anchor="end" font-family="Arial">{tick:.2f}%</text>')

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
