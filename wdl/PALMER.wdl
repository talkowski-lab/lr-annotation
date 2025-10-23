version 1.0

import "general/Structs.wdl"

workflow PALMER {
	input {
		File bam
		File bai

		String prefix
		String mode
		String MEI_type
		Array[String] contigs

		File ref_fa
		File ref_fai

		String utils_docker
		String palmer_docker
		RuntimeAttr? runtime_attr_split_bam
		RuntimeAttr? runtime_attr_run_palmer
		RuntimeAttr? runtime_attr_merge_palmer_outputs
	}

	call SplitBam {
		input:
			bam = bam,
			bai = bai,
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
				ref_fa = ref_fa,
				prefix = prefix,
				mode = mode,
				MEI_type = MEI_type,
				docker = palmer_docker,
				runtime_attr_override = runtime_attr_run_palmer
		}
	}

	call MergePALMEROutputs {
		input:
			calls_shards = RunPALMERShard.calls_shard,
			tsd_reads_shards = RunPALMERShard.tsd_reads_shard,
			prefix = prefix,
			docker = utils_docker,
			runtime_attr_override = runtime_attr_merge_palmer_outputs
	}

	output {
		File PALMER_calls = MergePALMEROutputs.calls
		File PALMER_tsd_reads = MergePALMEROutputs.tsd_reads
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
		max_retries: 1
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
		String MEI_type
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
			--type ~{MEI_type} \
			--mode ~{mode} \
			--output "~{prefix}" \
			--chr $chrom \
			--workdir "${dir}/${chrom}/"

		sed "s/$/\t~{MEI_type}/" ${chrom}/~{prefix}_calls.txt > ~{prefix}_calls_shard.txt || touch ~{prefix}_calls_shard.txt
		sed "s/$/\t~{MEI_type}/" ${chrom}/~{prefix}_tsd_reads.txt > ~{prefix}_tsd_reads_shard.txt || touch ~{prefix}_tsd_reads_shard.txt
	>>>

	output {
		File calls_shard = "~{prefix}_calls_shard.txt"
		File tsd_reads_shard = "~{prefix}_tsd_reads_shard.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 4,
		mem_gb: 12,
		disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) * 3 + 10,
		boot_disk_gb: 10,
		preemptible_tries: 1,
		max_retries: 1
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
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		head -n1 ~{calls_shards[0]} > ~{prefix}_calls.txt
		for f in ~{sep=' ' calls_shards}; do
			grep -v '^cluster_id' $f >> ~{prefix}_calls.txt || true
		done

		head -n1 ~{tsd_reads_shards[0]} > ~{prefix}_tsd_reads.txt
		for f in ~{sep=' ' tsd_reads_shards}; do
			grep -v '^cluster_id' $f >> ~{prefix}_tsd_reads.txt || true
		done
	>>>

	output {
		File calls = "~{prefix}_calls.txt"
		File tsd_reads = "~{prefix}_tsd_reads.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 2,
		disk_gb: ceil(size(calls_shards, "GB") + size(tsd_reads_shards, "GB")) * 2 + 10,
		boot_disk_gb: 10,
		preemptible_tries: 1,
		max_retries: 1
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
