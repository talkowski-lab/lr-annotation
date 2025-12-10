version 1.0

import "general/Structs.wdl"

workflow PALMER {
	input {
		File? bam_pat
		File? bai_pat
		Array[File]? override_palmer_calls_pat
		Array[File]? override_palmer_tsd_files_pat

		File? bam_mat
		File? bai_mat
		Array[File]? override_palmer_calls_mat
		Array[File]? override_palmer_tsd_files_mat

		String prefix
		String sample
		String mode
		Array[String] mei_types
		Array[String] contigs

		File ref_fa
		File ref_fai

		Array[String]? truvari_collapse_params

		String utils_docker
		String palmer_docker
		String annotate_palmer_docker

		RuntimeAttr? runtime_attr_split_bam
		RuntimeAttr? runtime_attr_run_palmer
		RuntimeAttr? runtime_attr_palmer_to_vcf
		RuntimeAttr? runtime_attr_truvari_collapse
		RuntimeAttr? runtime_attr_concat_sort_vcfs
	}

	scatter (idx in range(length(mei_types))) {
		String mei_type = mei_types[idx]
		String collapse_params = if defined(truvari_collapse_params) then select_first([truvari_collapse_params])[idx] else "--pctsize 0.9 --pctovl 0.9 --pctseq 0.9 --refdist 500"

		# Paternal haplotype
		if (!defined(override_palmer_calls_pat)) {
			call SplitBam as SplitBamPat {
				input:
					bam = select_first([bam_pat]),
					bai = select_first([bai_pat]),
					prefix = prefix + ".pat",
					contigs = contigs,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_split_bam
			}

			call RunPALMER as RunPALMERPat {
				input:
					bams = SplitBamPat.bams,
					bais = SplitBamPat.bais,
					prefix = prefix + ".pat",
					mode = mode,
					mei_type = mei_type,
					contigs = contigs,
					ref_fa = ref_fa,
					docker = palmer_docker,
					runtime_attr_override = runtime_attr_run_palmer
			}
		}

		File calls_file_pat = if defined (override_palmer_calls_pat) then select_first([override_palmer_calls_pat])[idx] else select_first([RunPALMERPat.calls])
		File tsd_file_pat = if defined (override_palmer_tsd_files_pat) then select_first([override_palmer_tsd_files_pat])[idx] else select_first([RunPALMERPat.tsd_reads])

		call ConvertPALMERToVcf as ConvertPALMERToVcfPat {
			input:
				PALMER_calls = calls_file_pat,
				PALMER_tsd_reads = tsd_file_pat,
				mei_type = mei_type,
				sample = sample,
				ref_fa = ref_fa,
				ref_fai = ref_fai,
				haplotype = "1|0",
				docker = annotate_palmer_docker,
				runtime_attr_override = runtime_attr_palmer_to_vcf
		}

		# Maternal haplotype
		if (!defined(override_palmer_calls_mat)) {
			call SplitBam as SplitBamMat {
				input:
					bam = select_first([bam_mat]),
					bai = select_first([bai_mat]),
					prefix = prefix + ".mat",
					contigs = contigs,
					docker = utils_docker,
					runtime_attr_override = runtime_attr_split_bam
			}

			call RunPALMER as RunPALMERMat {
				input:
					bams = SplitBamMat.bams,
					bais = SplitBamMat.bais,
					prefix = prefix + ".mat",
					mode = mode,
					mei_type = mei_type,
					contigs = contigs,
					ref_fa = ref_fa,
					docker = palmer_docker,
					runtime_attr_override = runtime_attr_run_palmer
			}
		}

		File calls_file_mat = if defined (override_palmer_calls_mat) then select_first([override_palmer_calls_mat])[idx] else select_first([RunPALMERMat.calls])
		File tsd_file_mat = if defined (override_palmer_tsd_files_mat) then select_first([override_palmer_tsd_files_mat])[idx] else select_first([RunPALMERMat.tsd_reads])

		call ConvertPALMERToVcf as ConvertPALMERToVcfMat {
			input:
				PALMER_calls = calls_file_mat,
				PALMER_tsd_reads = tsd_file_mat,
				mei_type = mei_type,
				sample = sample,
				ref_fa = ref_fa,
				ref_fai = ref_fai,
				haplotype = "0|1",
				docker = annotate_palmer_docker,
				runtime_attr_override = runtime_attr_palmer_to_vcf
		}

		# Merge haplotypes
		call TruvariCollapse {
			input:
				vcf_pat = ConvertPALMERToVcfPat.vcf,
				vcf_pat_idx = ConvertPALMERToVcfPat.vcf_idx,
				vcf_mat = ConvertPALMERToVcfMat.vcf,
				vcf_mat_idx = ConvertPALMERToVcfMat.vcf_idx,
				sample = sample,
				mei_type = mei_type,
				collapse_params = collapse_params,
				ref_fa = ref_fa,
				ref_fai = ref_fai,
				docker = utils_docker,
				runtime_attr_override = runtime_attr_truvari_collapse
		}
	}

	call ConcatSortVcfs {
		input:
			vcfs = TruvariCollapse.diploid_vcf,
			vcf_idxs = TruvariCollapse.diploid_vcf_idx,
			prefix = prefix,
			docker = utils_docker,
			runtime_attr_override = runtime_attr_concat_sort_vcfs
	}

	output {
		Array[File] PALMER_calls_pat = calls_file_pat
		Array[File] PALMER_tsd_reads_pat = tsd_file_pat
		Array[File] PALMER_calls_mat = calls_file_mat
		Array[File] PALMER_tsd_reads_mat = tsd_file_mat
		Array[File] PALMER_vcfs_pat = ConvertPALMERToVcfPat.vcf
		Array[File] PALMER_vcf_idxs_pat = ConvertPALMERToVcfPat.vcf_idx
		Array[File] PALMER_vcfs_mat = ConvertPALMERToVcfMat.vcf
		Array[File] PALMER_vcf_idxs_mat = ConvertPALMERToVcfMat.vcf_idx
		Array[File] PALMER_diploid_vcfs = TruvariCollapse.diploid_vcf
		Array[File] PALMER_diploid_vcf_idxs = TruvariCollapse.diploid_vcf_idx
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
		disk_gb: ceil(size(bam, "GB")) * 3 + 10,
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

task RunPALMER {
	input {
		Array[File] bams
		Array[File] bais
		String prefix
		String mode
		String mei_type
		Array[String] contigs
		File ref_fa
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -x

		dir=$(pwd)

		for bai in ~{sep=' ' bais}; do
			mv ${bai} ./
		done

		for bam_file in ~{sep=' ' bams}; do
			mv ${bam_file} ./
			bam=$(basename ${bam_file})
			chrom=$(echo $bam | sed 's/\.bam$//' | rev | cut -f1 -d '_' | rev)

			mkdir -p "${chrom}"
			/PALMER/PALMER \
				--input ${bam} \
				--ref_fa ~{ref_fa} \
				--ref_ver GRCh38 \
				--type ~{mei_type} \
				--mode ~{mode} \
				--output "~{prefix}" \
				--chr $chrom \
				--workdir "${dir}/${chrom}/" &
		done

		wait

		for chrom in ~{sep=' ' contigs}; do
			touch ${chrom}/~{prefix}_calls.txt
			touch ${chrom}/~{prefix}_TSD_reads.txt

			sed "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_calls.txt >> calls.txt
			sed "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_TSD_reads.txt >> TSD_reads.txt
		done

		head -n1 calls.txt > ~{prefix}_~{mei_type}_calls.txt
		grep -v '^cluster_id' calls.txt >> ~{prefix}_~{mei_type}_calls.txt || true

		head -n1 TSD_reads.txt > ~{prefix}_~{mei_type}_tsd_reads.txt
		grep -v '^cluster_id' TSD_reads.txt >> ~{prefix}_~{mei_type}_tsd_reads.txt || true
	>>>

	output {
		File calls = "~{prefix}_~{mei_type}_calls.txt"
		File tsd_reads = "~{prefix}_~{mei_type}_tsd_reads.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: length(bams),
		mem_gb: 6,
		disk_gb: ceil(size(bams, "GiB") + size(ref_fa, "GiB")) * 6,
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
		File PALMER_tsd_reads
		String mei_type
		String sample
		File ref_fa
		File ref_fai
		String haplotype
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		python /opt/gnomad-lr/scripts/mei/PALMER_to_vcf.py \
			--palmer_calls ~{PALMER_calls} \
			--palmer_tsd_reads ~{PALMER_tsd_reads} \
			--mei_type ~{mei_type} \
			--sample ~{sample} \
			--ref_fa ~{ref_fa} \
			--ref_fai ~{ref_fai} \
			--haplotype "~{haplotype}" \
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
		mem_gb: 4,
		disk_gb: ceil(size(PALMER_calls, "GB") + size(PALMER_tsd_reads, "GB") + size(ref_fa, "GB")) * 5 + 10,
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

task TruvariCollapse {
	input {
		File vcf_pat
		File vcf_pat_idx
		File vcf_mat
		File vcf_mat_idx
		String sample
		String mei_type
		String collapse_params
		File ref_fa
		File ref_fai
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		bcftools concat \
			-a \
			~{vcf_pat} \
			~{vcf_mat} \
		| bcftools sort \
			-Oz -o combined.vcf.gz
		
		tabix combined.vcf.gz

		truvari collapse \
			--reference ~{ref_fa} \
			-i combined.vcf.gz \
			-o ~{sample}.~{mei_type}.merged.vcf.gz \
			-c ~{sample}.~{mei_type}.collapsed.vcf.gz \
			--hap \
			~{collapse_params}

		bcftools sort \
			~{sample}.~{mei_type}.merged.vcf.gz \
			-Oz -o ~{sample}.~{mei_type}.merged.sorted.vcf.gz
		
		tabix ~{sample}.~{mei_type}.merged.sorted.vcf.gz
	>>>

	output {
		File diploid_vcf = "~{sample}.~{mei_type}.merged.sorted.vcf.gz"
		File diploid_vcf_idx = "~{sample}.~{mei_type}.merged.sorted.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: ceil(size(vcf_pat, "GB") + size(vcf_mat, "GB")) * 3 + 10,
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
		disk_gb: ceil(size(vcfs, "GB")) * 3 + 5,
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
