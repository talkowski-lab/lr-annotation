version 1.0

import "general/Structs.wdl"

workflow PALMER {
	input {
		File? bam
		File? bai

		String prefix
		String sample
		String? mode
		Array[String] mei_types
		Array[String]? contigs

		File? ref_fa
		File ref_fai

		Array[File]? override_palmer_calls
		Array[File]? override_palmer_tsd_files

		String utils_docker
		String palmer_docker

		RuntimeAttr? runtime_attr_split_bam
		RuntimeAttr? runtime_attr_run_palmer
		RuntimeAttr? runtime_attr_merge_palmer_outputs
		RuntimeAttr? runtime_attr_palmer_to_vcf
		RuntimeAttr? runtime_attr_concat_sort_vcfs
	}

	Boolean run_palmer = !defined(override_palmer_calls)

	if (run_palmer) {
		call SplitBam {
			input:
				bam = select_first([bam]),
				bai = select_first([bai]),
				prefix = prefix,
				contigs = select_first([contigs]),
				docker = utils_docker,
				runtime_attr_override = runtime_attr_split_bam
		}
	}

	scatter (idx in range(length(mei_types))) {
		String mei_type = mei_types[idx]

		if (run_palmer) {
			scatter (i in range(length(select_first([SplitBam.bams])))) {
				call RunPALMERShard {
					input:
						bam = select_first([SplitBam.bams])[i],
						bai = select_first([SplitBam.bais])[i],
						prefix = prefix,
						mode = select_first([mode]),
						mei_type = mei_type,
						ref_fa = select_first([ref_fa]),
						docker = palmer_docker,
						runtime_attr_override = runtime_attr_run_palmer
				}
			}

			call MergePALMEROutputs {
				input:
					calls_shards = select_first([RunPALMERShard.calls_shard]),
					tsd_reads_shards = select_first([RunPALMERShard.tsd_reads_shard]),
					prefix = prefix,
					mei_type = mei_type,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_merge_palmer_outputs
			}
		}

		File calls_file = select_first([MergePALMEROutputs.calls, select_first([override_palmer_calls])[idx]])
		File tsd_file = select_first([MergePALMEROutputs.tsd_reads, select_first([override_palmer_tsd_files])[idx]])

		call ConvertPALMERToVcf {
			input:
				PALMER_calls = calls_file,
				mei_type = mei_type,
				sample = sample,
				ref_fai = ref_fai,
				docker = palmer_docker,
				runtime_attr_override = runtime_attr_palmer_to_vcf
		}
	}

	call ConcatSortVcfs {
		input:
			vcfs = ConvertPALMERToVcf.vcf,
			vcf_idxs = ConvertPALMERToVcf.vcf_idx,
			prefix = prefix,
			docker = palmer_docker,
			runtime_attr_override = runtime_attr_concat_sort_vcfs
	}

	output {
		Array[File] PALMER_calls = calls_file
		Array[File] PALMER_tsd_reads = tsd_file
		Array[File] PALMER_vcfs = ConvertPALMERToVcf.vcf
		Array[File] PALMER_vcf_idxs = ConvertPALMERToVcf.vcf_idx
		File PALMER_combined_vcf = ConcatSortVcfs.vcf
		File PALMER_combined_vcf_idx = ConcatSortVcfs.vcf_idx
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
		disk_gb: ceil(size(bam, "GB") * 3) + 10,
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

		sed "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_calls.txt > ~{prefix}_calls_shard.txt || touch ~{prefix}_calls_shard.txt
		sed "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_TSD_reads.txt > ~{prefix}_tsd_reads_shard.txt || touch ~{prefix}_tsd_reads_shard.txt
	>>>

	output {
		File calls_shard = "~{prefix}_calls_shard.txt"
		File tsd_reads_shard = "~{prefix}_tsd_reads_shard.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 2,
		mem_gb: 4,
		disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) * 2 + 5,
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
		disk_gb: ceil(size(calls_shards, "GB") + size(tsd_reads_shards, "GB")) * 2 + 10,
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
		File PALMER_calls
		String mei_type
		String sample
		File ref_fai
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python /opt/gnomad-lr/scripts/mei/PALMER_to_vcf.py \
			--palmer_calls ~{PALMER_calls} \
			--mei_type ~{mei_type} \
			--sample ~{sample} \
			--ref_fai ~{ref_fai} \
			| bcftools sort -Oz \
			> ~{sample}.PALMER_calls.~{mei_type}.vcf.gz
		
		tabix ~{sample}.PALMER_calls.~{mei_type}.vcf.gz
	>>>

	output {
		File vcf = "~{sample}.PALMER_calls.~{mei_type}.vcf.gz"
		File vcf_idx = "~{sample}.PALMER_calls.~{mei_type}.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: 10*ceil(size(PALMER_calls, "GB") + size(ref_fai, "GB")) + 20,
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
		disk_gb: 5*ceil(size(vcfs, "GB")) + 5,
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
