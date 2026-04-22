version 1.0

import "../utils/Helpers.wdl"

workflow AnnotatePostProcess {
	input {
		File vcf
		File vcf_idx
		Array[String] contigs
		String prefix

		File ped

		Boolean filter_singletons = false

		String utils_docker

		RuntimeAttr? runtime_attr_subset
		RuntimeAttr? runtime_attr_post_process
		RuntimeAttr? runtime_attr_concat
	}

	Boolean single_contig = length(contigs) == 1

	scatter (contig in contigs) {
		if (!single_contig) {
			call Helpers.SubsetVcfToContig {
				input:
					vcf = vcf,
					vcf_idx = vcf_idx,
					contig = contig,
					prefix = prefix + "." + contig,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset
			}
		}

		File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
		File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

		call PostProcessVcf {
			input:
				vcf = contig_vcf,
				vcf_idx = contig_vcf_idx,
				contig = contig,
				ped = ped,
				filter_singletons = filter_singletons,
				prefix = prefix + "." + contig + ".post_process",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_post_process
		}
	}

	if (!single_contig) {
		call Helpers.ConcatVcfs {
			input:
				vcfs = PostProcessVcf.post_processed_vcf,
				vcf_idxs = PostProcessVcf.post_processed_vcf_idx,
				allow_overlaps = false,
				naive = true,
				prefix = prefix + ".post_processed",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_concat
		}
	}

	output {
		File post_processed_vcf = select_first([ConcatVcfs.concat_vcf, PostProcessVcf.post_processed_vcf[0]])
		File post_processed_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, PostProcessVcf.post_processed_vcf_idx[0]])
	}
}

task PostProcessVcf {
	input {
		File vcf
		File vcf_idx
		String contig
		File ped
		Boolean filter_singletons
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<CODE
import json
import math
import sys
import pysam
from math import comb

filter_singletons = ~{true="True" false="False" filter_singletons}
target_contig = "~{contig}"

sex_by_sample = {}
with open("~{ped}", "r") as handle:
	for line_number, line in enumerate(handle, start=1):
		fields = line.rstrip("\n").split("\t")
		if len(fields) < 5:
			print(f"PED line {line_number} has fewer than 5 columns", file=sys.stderr)
			sys.exit(1)

		sample_id = fields[1]
		sex_code = fields[4]
		if sample_id in sex_by_sample:
			print(f"Duplicate sample_id '{sample_id}' in PED file", file=sys.stderr)
			sys.exit(1)

		if sex_code == "1":
			sex_by_sample[sample_id] = "M"
		elif sex_code == "2":
			sex_by_sample[sample_id] = "F"
		else:
			sex_by_sample[sample_id] = None


def get_scalar(value):
	if isinstance(value, (list, tuple)):
		for item in value:
			if item is not None:
				return item
		return None
	return value


def genotype_field_length(n_alleles):
	return comb(n_alleles + 1, 2)


def is_missing(gt):
	return gt is None or all(allele is None for allele in gt)


def calculate_pl(ref_reads, alt_reads):
	total_reads = ref_reads + alt_reads
	if total_reads == 0:
		return (0, 0, 0)

	means = [0.03, 0.50, 0.97]
	priors = [0.33, 0.34, 0.33]

	ll_raw = []
	for prior, mean in zip(priors, means):
		ll = math.log(prior) + alt_reads * math.log(mean) + ref_reads * math.log(1.0 - mean)
		ll_raw.append(ll / math.log(10))

	max_ll = max(ll_raw)
	return tuple(int(round(-10 * (ll - max_ll))) for ll in ll_raw)


def calculate_gq(pls):
	return min(sorted(pls)[1], 99)


def prune_meis(record):
	allele_type = get_scalar(record.info.get("allele_type"))
	allele_length = get_scalar(record.info.get("allele_length"))
	if allele_type is None or allele_length is None:
		return

	length = abs(int(allele_length))
	if allele_type in {"alu_ins", "alu_del"} and (length < 250 or length > 350):
		record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
		if "SUB_FAMILY" in record.info:
			del record.info["SUB_FAMILY"]
	elif allele_type in {"sva_ins", "sva_del"} and (length < 1000 or length > 4000):
		record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
		if "SUB_FAMILY" in record.info:
			del record.info["SUB_FAMILY"]


def clear_genotype_fields(sample_data, n_alleles):
	gt = sample_data.get("GT")
	gt_length = len(gt) if gt is not None and len(gt) > 0 else 2
	sample_data["GT"] = tuple(None for _ in range(gt_length))
	if "DP" in sample_data:
		sample_data["DP"] = None
	if "GQ" in sample_data:
		sample_data["GQ"] = None
	if "AD" in sample_data:
		ad = sample_data.get("AD")
		ad_length = len(ad) if ad is not None and len(ad) > 0 else n_alleles
		sample_data["AD"] = tuple(None for _ in range(ad_length))
	if "PL" in sample_data:
		pl = sample_data.get("PL")
		pl_length = len(pl) if pl is not None and len(pl) > 0 else genotype_field_length(n_alleles)
		sample_data["PL"] = tuple(None for _ in range(pl_length))


def collapse_male_sex_chrom_gt(gt):
	if is_missing(gt):
		return gt

	left = gt[0] if len(gt) > 0 else None
	right = gt[1] if len(gt) > 1 else None

	if left is not None and left > 0 and not (right is not None and right > 0):
		return (left, None)
	if right is not None and right > 0 and not (left is not None and left > 0):
		return (None, right)
	if left is not None and left > 0:
		return (None, left)
	if right is not None and right > 0:
		return (None, right)
	return (0, None)


def update_haploid_likelihoods(sample_data, n_alleles):
	if "PL" not in sample_data and "GQ" not in sample_data:
		return


	ad = sample_data.get("AD")
	if ad is None or len(ad) != 2 or any(value is None for value in ad[:2]):
		if "PL" in sample_data:
			sample_data["PL"] = tuple(None for _ in range(genotype_field_length(n_alleles)))
		if "GQ" in sample_data:
			sample_data["GQ"] = None
		return

	pls = calculate_pl(ad[0], ad[1])
	if "PL" in sample_data:
		sample_data["PL"] = pls
	if "GQ" in sample_data:
		sample_data["GQ"] = calculate_gq(pls)


def shortest_motif_length(record):
	motifs = record.info.get("MOTIFS")
	if motifs is None:
		return None

	motif_values = []
	if isinstance(motifs, (list, tuple)):
		raw_values = motifs
	else:
		raw_values = [motifs]

	for value in raw_values:
		if value is None:
			continue
		motif_values.extend(part for part in str(value).split(",") if part and part != ".")

	if not motif_values:
		return None
	return min(len(motif) for motif in motif_values)


def has_single_read_support(record):
	if not filter_singletons:
		return False
    
    if 'AC' in record.info and len(record.alts) == 1 and record.info['AC'][0] > 2:
        return False

	alt_depths = []
	for sample_data in record.samples.values():
		ad = sample_data.get("AD")
		if ad is None or len(ad) < 2:
			continue

		alt_depth = 0
		has_alt_depth = False
		for value in ad[1:]:
			if value is not None:
				alt_depth += value
				has_alt_depth = True

		if has_alt_depth and alt_depth > 0:
			alt_depths.append(alt_depth)

	return len(alt_depths) == 1 and alt_depths[0] == 1

# Update header
vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header.copy()
if "HOMOPOLYMER_TRV" not in header.filters:
	header.filters.add("HOMOPOLYMER_TRV", None, None, "TRV call where the shortest motif has length 1")
if filter_singletons and "SINGLE_READ_SUPPORT" not in header.filters:
	header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by only one read in a single sample")

# Set up variables
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
header_info = vcf_out.header.info
samples_in_vcf = list(vcf_in.header.samples)
iterator = vcf_in.fetch(target_contig)

# Iterate through records
for record in iterator:
	record.translate(vcf_out.header)
	prune_meis(record)

	# Clear genotype fields for females on chrY
	n_alleles = len(record.alleles)
	if record.chrom == "chrY":
		for sample in samples_in_vcf:
			if sex_by_sample.get(sample) == "F":
				clear_genotype_fields(record.samples[sample], n_alleles)

	# Set GTs to hemizygous for males on chrX/chrY
	if record.chrom in {"chrX", "chrY"}:
		for sample in samples_in_vcf:
			if sex_by_sample.get(sample) != "M":
				continue

			sample_data = record.samples[sample]
			gt = sample_data.get("GT")
			if is_missing(gt):
				continue

			sample_data["GT"] = collapse_male_sex_chrom_gt(gt)
			update_haploid_likelihoods(sample_data, n_alleles)

	# Flag homopolymer TR variants
	if get_scalar(record.info.get("allele_type")) == "trv" and shortest_motif_length(record) == 1:
		record.filter.add("HOMOPOLYMER_TRV")

	# Flag variants with single read support
	if has_single_read_support(record):
		record.filter.add("SINGLE_READ_SUPPORT")

	vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

		tabix -p vcf -f ~{prefix}.vcf.gz
	>>>

	output {
		File post_processed_vcf = "~{prefix}.vcf.gz"
		File post_processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: ceil(size(vcf, "GB")) + 5,
		disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
		boot_disk_gb: 10,
		preemptible_tries: 1,
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
