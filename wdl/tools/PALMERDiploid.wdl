version 1.0

import "../utils/Structs.wdl"

workflow PALMERDiploid {
	input {
		File? bam
		File? bai
		Array[File]? override_palmer_calls
		Array[File]? override_palmer_tsd_files

		String prefix
		String sample
		String mode
		Array[String] mei_types
		Array[String] contigs

		File ref_fa
		File ref_fai

		String utils_docker
		String palmer_docker
		String annotate_palmer_docker

		RuntimeAttr? runtime_attr_split_bam
		RuntimeAttr? runtime_attr_run_palmer
		RuntimeAttr? runtime_attr_merge_palmer_outputs
		RuntimeAttr? runtime_attr_palmer_to_vcf
		RuntimeAttr? runtime_attr_concat_sort_vcfs
	}

	scatter (idx in range(length(mei_types))) {
		String mei_type = mei_types[idx]

		if (!defined(override_palmer_calls)) {
			call SplitBam {
				input:
					bam = select_first([bam]),
					bai = select_first([bai]),
					prefix = prefix,
					contigs = contigs,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_split_bam
			}

			scatter (i in range(length(SplitBam.bams))) {
				call RunPALMERShard {
					input:
						bam = SplitBam.bams[i],
						bai = SplitBam.bais[i],
						prefix = prefix,
						mode = mode,
						mei_type = mei_type,
						ref_fa = ref_fa,
						docker = palmer_docker,
						runtime_attr_override = runtime_attr_run_palmer
				}
			}

			call MergePALMEROutputs {
				input:
					calls_shards = RunPALMERShard.calls_shard,
					tsd_reads_shards = RunPALMERShard.tsd_reads_shard,
					prefix = prefix,
					mei_type = mei_type,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_merge_palmer_outputs
			}
		}

		if (defined(override_palmer_calls)) {
			call AddMeiTypeColumn as AddMeiTypeToCallsOverride {
				input:
					input_file = select_first([override_palmer_calls])[idx],
					mei_type = mei_type,
					prefix = prefix,
					file_type = "calls",
					docker = utils_docker
			}
		}

		if (defined(override_palmer_tsd_files)) {
			call AddMeiTypeColumn as AddMeiTypeToTsdOverride {
				input:
					input_file = select_first([override_palmer_tsd_files])[idx],
					mei_type = mei_type,
					prefix = prefix,
					file_type = "tsd_reads",
					docker = utils_docker
			}
		}

		File calls_file = if defined (override_palmer_calls) then select_first([AddMeiTypeToCallsOverride.output_file]) else select_first([MergePALMEROutputs.calls])
		File tsd_file = if defined (override_palmer_tsd_files) then select_first([AddMeiTypeToTsdOverride.output_file]) else select_first([MergePALMEROutputs.tsd_reads])

		call ConvertPALMERToVcf {
			input:
				palmer_calls = calls_file,
				palmer_tsd_reads = tsd_file,
				mei_type = mei_type,
				sample = sample,
				ref_fa = ref_fa,
				ref_fai = ref_fai,
				docker = annotate_palmer_docker,
				runtime_attr_override = runtime_attr_palmer_to_vcf
		}
	}

	call ConcatSortVcfs {
		input:
			vcfs = ConvertPALMERToVcf.vcf,
			vcf_idxs = ConvertPALMERToVcf.vcf_idx,
			prefix = prefix,
			docker = utils_docker,
			runtime_attr_override = runtime_attr_concat_sort_vcfs
	}

	output {
		Array[File] palmer_calls = calls_file
		Array[File] palmer_tsd_reads = tsd_file
		Array[File] palmer_diploid_vcfs = ConvertPALMERToVcf.vcf
		Array[File] palmer_diploid_vcfs_idxs = ConvertPALMERToVcf.vcf_idx

		File palmer_combined_vcf = ConcatSortVcfs.vcf
		File palmer_combined_vcf_idx = ConcatSortVcfs.vcf_idx
	}
}

task SplitBam {
	input {
		File bam
		File bai
		String prefix
		Array[String] contigs
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		for contig in ~{sep=" " contigs}
		do
			samtools view -b ~{bam} $contig > ~{prefix}_${contig}.bam
			samtools index ~{prefix}_${contig}.bam
		done
	>>>

	output {
		Array[File] bams = glob("~{prefix}_*bam")
		Array[File] bais = glob("~{prefix}_*bai")
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 3 * ceil(size(bam, "GB")) + 10,
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

task RunPALMERShard {
	input {
		File bam
		File bai
		String prefix
		String mode
		String mei_type
		File ref_fa
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		dir=$(pwd)

		mv ~{bam} ./
		mv ~{bai} ./
		bam_base=$(basename ~{bam})
		chrom=$(echo $bam_base | sed 's/\.bam$//' | rev | cut -f1 -d '_' | rev)

		mkdir -p "${chrom}"
		/PALMER/PALMER \
			--input ${bam_base} \
			--ref_fa ~{ref_fa} \
			--ref_ver GRCh38 \
			--type ~{mei_type} \
			--mode ~{mode} \
			--output "~{prefix}" \
			--chr $chrom \
			--workdir "${dir}/${chrom}/"

		sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_calls.txt
		sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_TSD_reads.txt
		mv ${chrom}/~{prefix}_calls.txt ~{prefix}_calls_shard.txt
		mv ${chrom}/~{prefix}_TSD_reads.txt ~{prefix}_tsd_reads_shard.txt
	>>>

	output {
		File calls_shard = "~{prefix}_calls_shard.txt"
		File tsd_reads_shard = "~{prefix}_tsd_reads_shard.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 4,
		boot_disk_gb: 5,
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

task MergePALMEROutputs {
	input {
		Array[File] calls_shards
		Array[File] tsd_reads_shards
		String prefix
		String mei_type
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		head -n1 ~{calls_shards[0]} > ~{prefix}_~{mei_type}_calls.txt
		for f in ~{sep=' ' calls_shards}; do
			grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_calls.txt || true
		done

		head -n1 ~{tsd_reads_shards[0]} > ~{prefix}_~{mei_type}_tsd_reads.txt
		for f in ~{sep=' ' tsd_reads_shards}; do
			grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_tsd_reads.txt || true
		done
	>>>

	output {
		File calls = "~{prefix}_~{mei_type}_calls.txt"
		File tsd_reads = "~{prefix}_~{mei_type}_tsd_reads.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: 2 * ceil(size(calls_shards, "GB") + size(tsd_reads_shards, "GB")) + 10,
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

task ConvertPALMERToVcf {
	input {
		File palmer_calls
		File palmer_tsd_reads
		String mei_type
		String sample
		File ref_fa
		File ref_fai
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python /opt/gnomad-lr/scripts/mei/PALMER_to_vcf.py \
			--palmer_calls ~{palmer_calls} \
			--palmer_tsd_reads ~{palmer_tsd_reads} \
			--mei_type ~{mei_type} \
			--sample ~{sample} \
			--ref_fa ~{ref_fa} \
			--ref_fai ~{ref_fai} \
			--haplotype "1/1" \
			| bcftools sort -Oz \
			> ~{sample}.palmer_calls.~{mei_type}.vcf.gz
		
		tabix ~{sample}.palmer_calls.~{mei_type}.vcf.gz
	>>>

	output {
		File vcf = "~{sample}.palmer_calls.~{mei_type}.vcf.gz"
		File vcf_idx = "~{sample}.palmer_calls.~{mei_type}.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: 5 * ceil(size(palmer_calls, "GB") + size(palmer_tsd_reads, "GB") + size(ref_fa, "GB")) + 10,
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

task AddMeiTypeColumn {
	input {
		File input_file
		String mei_type
		String prefix
		String file_type
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		sed "s/$/\t~{mei_type}/" ~{input_file} > ~{prefix}_~{mei_type}_~{file_type}.txt
	>>>

	output {
		File output_file = "~{prefix}_~{mei_type}_~{file_type}.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: 2 * ceil(size(input_file, "GB")) + 10,
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

task ConcatSortVcfs {
	input {
		Array[File] vcfs
		Array[File] vcf_idxs
		String prefix
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		bcftools concat -a ~{sep=" " vcfs} | bcftools sort -Oz > ~{prefix}.vcf.gz
		tabix ~{prefix}.vcf.gz
	>>>

	output {
		File vcf = "~{prefix}.vcf.gz"
		File vcf_idx = "~{prefix}.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: 3 * ceil(size(vcfs, "GB")) + 5,
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
