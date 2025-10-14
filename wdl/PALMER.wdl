version 1.0

import "Structs.wdl"

workflow PALMER {
	input {
		File bam
		File bai
		File ref_fa
		File ref_fai
		String prefix
		String mode
		String MEI_type
		Array[String] contigs

		String utils_docker
		String palmer_docker
	}

	call split_bam {
		input:
			bam = bam,
			bai = bai,
			prefix = prefix,
			contigs = contigs,
			docker = utils_docker
	}

	call PALMER_WGS {
		input:
			bams = split_bam.bams,
			bais = split_bam.bais,
			ref_fa = ref_fa,
			prefix = prefix,
			mode = mode,
			MEI_type = MEI_type,
			contigs = contigs,
			docker = palmer_docker
	}

	output {
		File PALMER_calls = PALMER_WGS.calls
		File PALMER_TSD_reads = PALMER_WGS.TSD_reads
	}
}

task split_bam {
	input {
		File bam
		File bai
		String prefix
		Array[String] contigs
		
		String docker
		RuntimeAttr? runtime_attr_override
	}

	Int disk_size = ceil(size(bam, "GB")) * 3

	command <<<
		set -euxo pipefail

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
		cpu_cores:          1,
		mem_gb:             4,
		disk_gb:            disk_size,
		boot_disk_gb:       10,
		preemptible_tries:  2,
		max_retries:        1,
		docker:             docker
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
		memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
		disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
		preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
		docker:                 select_first([runtime_attr.docker,            default_attr.docker])
	}
}

task PALMER_WGS {
	input {
		Array[File] bams
		Array[File] bais
		File ref_fa
		String prefix
		String mode
		String MEI_type
		Array[String] contigs

		String docker
		RuntimeAttr? runtime_attr_override
	}

	parameter_meta {
		bams:             "aligned bams, can contain long reads or assembled contigs"
		bais:             "indices accompanying the BAM"
		ref_fa:           "fasta of reference used"
		prefix:           "output file prefix"
		mode:             "either raw or asm"
		MEI_type:         "type of MEI to call, options are LINE, ALU, SVA, HERVK"
	}

	Int disk_size = ceil(size(bams, "GB") + size(ref_fa, "GB")) * 6
	Int cpu_num = length(bams)

	command <<<
		set -x

		dir=$(pwd)

		for bai in ~{sep=' ' bais}; do
			mv ${bai} ./
		done

		chroms=""

		for bam_file in ~{sep=' ' bams}; do
			mv ${bam_file} ./
			bam=$(basename ${bam_file})

			# taking advantage of naming scheme defined in split_bam to get name of chrom
			chrom=$(echo $bam|sed 's/\.bam$//'|rev|cut -f1 -d '_'|rev)

			mkdir -p "${chrom}"
			/PALMER/PALMER --input ${bam} \
					--ref_fa ~{ref_fa} \
					--ref_ver GRCh38 \
					--type ~{MEI_type} \
					--mode ~{mode} \
					--output "~{prefix}" \
					--chr $chrom \
					--workdir "${dir}/${chrom}/" &
		done

		wait # wait for all background PALMER jobs to finish

		# concat results and annotate with MEI type
		for chrom in ~{sep=' ' contigs}; do
			touch ${chrom}/~{prefix}_calls.txt
			touch ${chrom}/~{prefix}_TSD_reads.txt

			sed "s/$/\t~{MEI_type}/" ${chrom}/~{prefix}_calls.txt >> calls.txt
			sed "s/$/\t~{MEI_type}/" ${chrom}/~{prefix}_TSD_reads.txt >> TSD_reads.txt
		done

		head -n1 calls.txt > ~{prefix}_calls.txt
		grep -v '^cluster_id' calls.txt >> ~{prefix}_calls.txt

		head -n1 TSD_reads.txt > ~{prefix}_TSD_reads.txt
		grep -v '^cluster_id' TSD_reads.txt >> ~{prefix}_TSD_reads.txt
	>>>

	output {
		File calls = "~{prefix}_calls.txt"
		File TSD_reads = "~{prefix}_TSD_reads.txt"
	}

	RuntimeAttr default_attr = object {
		cpu_cores:          cpu_num,
		mem_gb:             4,         # PALMER doesn't use a lot of memory
		disk_gb:            disk_size,
		boot_disk_gb:       10,
		preemptible_tries:  3,
		max_retries:        0,
		docker:             docker
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
		memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
		disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
		preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
		docker:                 select_first([runtime_attr.docker,            default_attr.docker])
	}
}
