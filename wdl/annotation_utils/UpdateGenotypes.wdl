version 1.0

import "../utils/Helpers.wdl"

workflow UpdateGenotypes {
	input {
		File base_vcf
		File base_vcf_idx
		File? genotyped_vcf
		File? genotyped_vcf_idx
		File ped

		Array[String] contigs
		String prefix

		Boolean transfer_genotypes = false

		String utils_docker

		RuntimeAttr? runtime_attr_subset_base
		RuntimeAttr? runtime_attr_subset_genotyped
		RuntimeAttr? runtime_attr_update_genotypes
		RuntimeAttr? runtime_attr_concat
	}

	Boolean single_contig = length(contigs) == 1

	scatter (contig in contigs) {
		if (!single_contig) {
			call Helpers.SubsetVcfToContig as SubsetBase {
				input:
					vcf = base_vcf,
					vcf_idx = base_vcf_idx,
					contig = contig,
					prefix = prefix + "." + contig + ".base",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_base
			}
		}

		if (transfer_genotypes && !single_contig) {
			call Helpers.SubsetVcfToContig as SubsetGenotyped {
				input:
					vcf = select_first([genotyped_vcf]),
					vcf_idx = select_first([genotyped_vcf_idx]),
					contig = contig,
					prefix = prefix + "." + contig + ".genotyped",
					docker = utils_docker,
					runtime_attr_override = runtime_attr_subset_genotyped
			}
		}

		File contig_base_vcf = select_first([SubsetBase.subset_vcf, base_vcf])
		File contig_base_vcf_idx = select_first([SubsetBase.subset_vcf_idx, base_vcf_idx])
		File transfer_source_vcf = if transfer_genotypes then select_first([SubsetGenotyped.subset_vcf, select_first([genotyped_vcf])]) else contig_base_vcf
		File transfer_source_vcf_idx = if transfer_genotypes then select_first([SubsetGenotyped.subset_vcf_idx, select_first([genotyped_vcf_idx])]) else contig_base_vcf_idx

		call UpdateContigGenotypes {
			input:
				base_vcf = contig_base_vcf,
				base_vcf_idx = contig_base_vcf_idx,
				genotyped_vcf = transfer_source_vcf,
				genotyped_vcf_idx = transfer_source_vcf_idx,
				ped = ped,
				contig = contig,
				transfer_genotypes = transfer_genotypes,
				prefix = prefix + "." + contig + ".updated",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_update_genotypes
		}
	}

	if (!single_contig) {
		call Helpers.ConcatVcfs {
			input:
				vcfs = UpdateContigGenotypes.updated_vcf,
				vcf_idxs = UpdateContigGenotypes.updated_vcf_idx,
				allow_overlaps = false,
				naive = true,
				prefix = prefix + ".updated",
				docker = utils_docker,
				runtime_attr_override = runtime_attr_concat
		}
	}

	output {
		File updated_vcf = select_first([ConcatVcfs.concat_vcf, UpdateContigGenotypes.updated_vcf[0]])
		File updated_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, UpdateContigGenotypes.updated_vcf_idx[0]])
	}
}

task UpdateContigGenotypes {
	input {
		File base_vcf
		File base_vcf_idx
		File genotyped_vcf
		File genotyped_vcf_idx
		File ped
		String contig

		Boolean transfer_genotypes

		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python3 <<CODE
import sys

import pysam


def parse_ped(path):
	sex_by_sample = {}
	with open(path, "r") as handle:
		for line_number, line in enumerate(handle, start=1):
			fields = line.strip().split()
			if not fields:
				continue

			sample_id = fields[1]
			sex_code = fields[4]
			if sex_code == "1":
				sex_by_sample[sample_id] = "M"
			elif sex_code == "2":
				sex_by_sample[sample_id] = "F"
			else:
				sex_by_sample[sample_id] = None
	return sex_by_sample


def is_heterozygous(gt):
	if gt is None:
		return False
	called = [allele for allele in gt if allele is not None]
	return len(called) >= 2 and len(set(called)) > 1


def clear_format_fields(sample_data):
	gt = sample_data.get("GT")
	gt_length = len(gt) if gt is not None and len(gt) > 0 else 2
	sample_data["GT"] = tuple(None for _ in range(gt_length))
	sample_data.phased = False

	for fmt_key in list(sample_data.keys()):
		if fmt_key == "GT":
			continue
		value = sample_data.get(fmt_key)
		if isinstance(value, tuple):
			sample_data[fmt_key] = tuple(None for _ in range(len(value)))
		else:
			sample_data[fmt_key] = None


def right_align_unphased(gt):
	if gt is None:
		return gt
	return tuple(sorted(gt, key=lambda allele: (allele is not None, allele if allele is not None else -1)))


def make_male_hemizygous(gt, phased):
	if gt is None:
		return gt

	alleles = list(gt)
	called_positions = [index for index, allele in enumerate(alleles) if allele is not None]
	if len(called_positions) <= 1:
		return tuple(alleles)

	alt_positions = [index for index, allele in enumerate(alleles) if allele is not None and allele > 0]

	if phased:
		if len(alt_positions) == 1:
			keep_index = alt_positions[0]
		elif alt_positions:
			keep_index = alt_positions[-1]
		else:
			keep_index = called_positions[-1]

		new_gt = [None] * len(alleles)
		new_gt[keep_index] = alleles[keep_index]
		return tuple(new_gt)

	if alt_positions:
		keep_allele = alleles[alt_positions[-1]]
	else:
		keep_allele = alleles[called_positions[-1]]

	new_gt = [None] * len(alleles)
	new_gt[-1] = keep_allele
	return right_align_unphased(tuple(new_gt))


def find_match(record, genotyped_reader):
	for candidate in genotyped_reader.fetch(record.chrom, record.start, record.stop):
		if candidate.id == record.id and candidate.ref == record.ref and candidate.alts == record.alts:
			return candidate
	return None


transfer_genotypes = ~{true="True" false="False" transfer_genotypes}
target_contig = "~{contig}"
sex_by_sample = parse_ped("~{ped}")

base_reader = pysam.VariantFile("~{base_vcf}")
genotyped_reader = pysam.VariantFile("~{genotyped_vcf}") if transfer_genotypes else None

output_writer = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=base_reader.header)
base_samples = list(base_reader.header.samples)
shared_samples = set(base_samples)
if transfer_genotypes:
	shared_samples &= set(genotyped_reader.header.samples)

try:
	iterator = base_reader.fetch(target_contig)
except ValueError:
	iterator = iter(())

for record in iterator:
	record.translate(output_writer.header)

	# Transfer GT field
	if transfer_genotypes:
		match = find_match(record, genotyped_reader)
		if match is not None:
			for sample in shared_samples:
				base_gt = record.samples[sample].get("GT")
				if is_heterozygous(base_gt):
					record.samples[sample]["GT"] = match.samples[sample].get("GT")
					record.samples[sample].phased = match.samples[sample].phased

	for sample in base_samples:
		sample_data = record.samples[sample]
		sample_sex = sex_by_sample.get(sample)

		# Clear format fields for females on chrY
		if record.chrom == "chrY" and sample_sex == "F":
			clear_format_fields(sample_data)
			continue

		# Make male calls hemizygous on chrX & chrY
		if record.chrom in {"chrX", "chrY"} and sample_sex == "M":
			sample_data["GT"] = make_male_hemizygous(sample_data.get("GT"), sample_data.phased)

		# Right align unphased calls
		if not sample_data.phased:
			sample_data["GT"] = right_align_unphased(sample_data.get("GT"))

	output_writer.write(record)

base_reader.close()
output_writer.close()
if genotyped_reader is not None:
	genotyped_reader.close()
CODE

		tabix -p vcf -f ~{prefix}.vcf.gz
	>>>

	output {
		File updated_vcf = "~{prefix}.vcf.gz"
		File updated_vcf_idx = "~{prefix}.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 5 * ceil(size(base_vcf, "GB") + size(genotyped_vcf, "GB")) + 25,
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
