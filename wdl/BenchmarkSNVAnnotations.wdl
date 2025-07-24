version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers
import "ParseVcfSharded.wdl" as ParseVcf

workflow BenchmarkSNVAnnotations {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth
        File vcf_truth_index
        File primary_contigs_list
        String pipeline_docker
        String prefix

        Int? variants_per_shard
        
        RuntimeAttr? runtime_attr_get_vep_format
        RuntimeAttr? runtime_attr_subset_eval_contig
        RuntimeAttr? runtime_attr_subset_truth_contig
        RuntimeAttr? runtime_attr_split_eval_vcf
        RuntimeAttr? runtime_attr_split_truth_vcf
        RuntimeAttr? runtime_attr_parse_eval_vcf
        RuntimeAttr? runtime_attr_parse_truth_vcf
        RuntimeAttr? runtime_attr_concat_eval_tsv
        RuntimeAttr? runtime_attr_concat_truth_tsv
        RuntimeAttr? runtime_attr_benchmark
        RuntimeAttr? runtime_attr_merge_summaries
        RuntimeAttr? runtime_attr_merge_tarballs
    }

    Array[String] contigs = read_lines(primary_contigs_list)
    Int variants_per_shard_eff = select_first([variants_per_shard, 100000000])

    call GetVepFormat as GetVepFormatEval {
        input:
            vcf = vcf_eval,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_get_vep_format
    }

    call GetVepFormat as GetVepFormatTruth {
        input:
            vcf = vcf_truth,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_get_vep_format
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetEvalVcf {
            input:
                vcf = vcf_eval,
                vcf_index = vcf_eval_index,
                contig = contig,
                prefix = "~{prefix}.eval",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_eval_contig
        }

        call Helpers.SubsetVcfToContig as SubsetTruthVcf {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_index,
                contig = contig,
                prefix = "~{prefix}.truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_truth_contig
        }

        call ParseVcf.ParseVcfSharded as ParseEvalVcf {
            input:
                vcf = SubsetEvalVcf.subset_vcf,
                vcf_index = SubsetEvalVcf.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.eval",
                variants_per_shard = variants_per_shard_eff,
                pipeline_docker = pipeline_docker,
                runtime_attr_split_vcf = runtime_attr_split_eval_vcf,
                runtime_attr_parse_vcf = runtime_attr_parse_eval_vcf,
                runtime_attr_concat_tsv = runtime_attr_concat_eval_tsv
        }

        call ParseVcf.ParseVcfSharded as ParseTruthVcf {
            input:
                vcf = SubsetTruthVcf.subset_vcf,
                vcf_index = SubsetTruthVcf.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.truth",
                variants_per_shard = variants_per_shard_eff,
                pipeline_docker = pipeline_docker,
                runtime_attr_split_vcf = runtime_attr_split_truth_vcf,
                runtime_attr_parse_vcf = runtime_attr_parse_truth_vcf,
                runtime_attr_concat_tsv = runtime_attr_concat_truth_tsv
        }

        call BenchmarkContig {
            input:
                eval_tsv = ParseEvalVcf.output_tsv,
                truth_tsv = ParseTruthVcf.output_tsv,
                eval_vep_format = GetVepFormatEval.vep_format,
                truth_vep_format = GetVepFormatTruth.vep_format,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_benchmark
        }
    }

    call MergeSummaries {
        input:
            summary_files = BenchmarkContig.summary_file,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_summaries
    }

    call MergePlotTarballs {
        input:
            tarballs = BenchmarkContig.plot_tarball,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_tarballs
    }

    output {
        File merged_summary = MergeSummaries.merged_file
        File merged_plot_tarball = MergePlotTarballs.merged_tarball
    }
}

task GetVepFormat {
    input {
        File vcf
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/get_vep_format.py ~{vcf} > vep_format.txt
    >>>

    output {
        String vep_format = read_string("vep_format.txt")
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task BenchmarkContig {
    input {
        File eval_tsv
        File truth_tsv
        String eval_vep_format
        String truth_vep_format
        String contig
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/benchmark_snv_annotations.py ~{eval_tsv} ~{truth_tsv} "~{eval_vep_format}" "~{truth_vep_format}" ~{contig} ~{prefix}
    >>>

    output {
        File summary_file = "~{prefix}_summary.txt"
        File plot_tarball = "~{prefix}.tar.gz"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(eval_tsv, "GB") + size(truth_tsv, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeSummaries {
    input {
        Array[File] summary_files
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        cat ~{sep=' ' summary_files} > ~{prefix}.summary.txt
    >>>

    output {
        File merged_file = "~{prefix}.summary.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(summary_files, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergePlotTarballs {
    input {
        Array[File] tarballs
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        mkdir -p plots/AF_plots
        mkdir -p plots/VEP_plots

        for tarball in ~{sep=' ' tarballs}; do
            tar -xzf $tarball --strip-components=1 -C plots/
        done

        tar -czf ~{prefix}.plots.tar.gz plots/

    >>>

    output {
        File merged_tarball = "~{prefix}.plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(tarballs, "GB") * 1.5) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
