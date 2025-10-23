version 1.0

import "general/Structs.wdl"

workflow PAV {
	input {
		Array[File] mat_haplotypes
		Array[File] pat_haplotypes
		Array[String] sample_ids

		File ref_fasta
		File ref_fasta_fai

		String pav_docker

		RuntimeAttr? runtime_attr_run_pav
	}

	call CallPAV {
		input:
			mat_haplotypes = mat_haplotypes,
			pat_haplotypes = pat_haplotypes,
			sample_ids = sample_ids,
			ref_fasta = ref_fasta,
			ref_fasta_fai = ref_fasta_fai,
			pav_docker = pav_docker,
			runtime_attr_override = runtime_attr_run_pav
	}

	output {
		File pav_results_tarball = CallPAV.results_tar
		File pav_log_tarball = CallPAV.log_tar
	}
}

task CallPAV {
	input {
		Array[File] mat_haplotypes
		Array[File] pat_haplotypes
		Array[String] sample_ids

		File ref_fasta
		File ref_fasta_fai

		String pav_docker

		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size(mat_haplotypes, "GiB") + size(pat_haplotypes, "GiB") + size(ref_fasta, "GiB")
	Int disk_size = ceil(input_size * 3) + 20

	RuntimeAttr default_attr = object {
		cpu_cores: 16,
		mem_gb: 64,
		disk_gb: disk_size,
		boot_disk_gb: 10,
		preemptible_tries: 1,
		max_retries: 1
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	Int effective_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

	command <<<
		set -euo pipefail

		ln -s ~{ref_fasta} ref.fa
		ln -s ~{ref_fasta_fai} ref.fa.fai

		python3 <<'CODE'
import os
import json

config = {"reference": "ref.fa"}
with open("config.json", "w") as f:
	json.dump(config, f)

os.makedirs("asms", exist_ok=True)

sample_ids = "~{sep=' ' sample_ids}".split(' ')
mat_files = "~{sep=' ' mat_haplotypes}".split(' ')
pat_files = "~{sep=' ' pat_haplotypes}".split(' ')

if len(sample_ids) != len(mat_files) or len(sample_ids) != len(pat_files):
	raise ValueError(f"Input array lengths must match: {len(sample_ids)} samples, {len(mat_files)} maternal, {len(pat_files)} paternal")

with open("assemblies.tsv", "w") as f:
	f.write("NAME\tHAP_mat\tHAP_pat\n")
	
	for i, sample_id in enumerate(sample_ids):
		mat_link_rel = f"asms/{sample_id}_mat.fa.gz"
		pat_link_rel = f"asms/{sample_id}_pat.fa.gz"

		os.symlink(mat_files[i], mat_link_rel)
		os.symlink(pat_files[i], pat_link_rel)

		mat_link_abs = os.path.abspath(mat_link_rel)
		pat_link_abs = os.path.abspath(pat_link_rel)
		
		f.write(f"{sample_id}\t{mat_link_abs}\t{pat_link_abs}\n")
CODE

		python3 -m pav3 call --cores ~{effective_cpu}

		tar -zcf pav_results.tar.gz results
		tar -zcf pav_log.tar.gz log
	>>>

	output {
		File results_tar = "pav_results.tar.gz"
		File log_tar = "pav_log.tar.gz"
	}

	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: pav_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}
