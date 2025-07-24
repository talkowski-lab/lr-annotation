version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers
import "BedtoolsClosestSV.wdl" as Bedtools

workflow BenchmarkAnnotations {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth
        File vcf_truth_index
        File vcf_sv_truth
        File vcf_sv_truth_index
        File ref_fasta
        File ref_fasta_index
        
        File primary_contigs_list
        String pipeline_docker
        String prefix
        Boolean create_benchmarks

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_benchmark
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_merge_summaries
        RuntimeAttr? runtime_attr_merge_tarballs
        RuntimeAttr? runtime_attr_bedtools_closest
    }

    Array[String] contigs = read_lines(primary_contigs_list)

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetEvalVcf {
            input:
                vcf = vcf_eval,
                vcf_index = vcf_eval_index,
                contig = contig,
                prefix = "~{prefix}.eval",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        call Helpers.SubsetVcfToContig as SubsetTruthVcf {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_index,
                contig = contig,
                prefix = "~{prefix}.truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        call Helpers.SubsetVcfToContig as SubsetSVTruthVcf {
            input:
                vcf = vcf_sv_truth,
                vcf_index = vcf_sv_truth_index,
                contig = contig,
                prefix = "~{prefix}.sv_truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        call Bedtools.BedtoolsClosestSV as BedtoolsClosest {
            input:
                vcf_eval = SubsetEvalVcf.subset_vcf,
                vcf_truth = SubsetSVTruthVcf.subset_vcf, # Using SV truth set
                prefix = "~{prefix}.~{contig}",
                sv_pipeline_docker = pipeline_docker,
                runtime_attr_split_vcf = runtime_attr_subset_contig,
                runtime_attr_bedtools_closest = runtime_attr_bedtools_closest
        }

        call BenchmarkContig {
            input:
                vcf_eval = SubsetEvalVcf.subset_vcf,
                vcf_truth = SubsetTruthVcf.subset_vcf,
                sv_truth = SubsetSVTruthVcf.subset_vcf,
                closest_bed = BedtoolsClosest.closest_bed,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                create_benchmarks = create_benchmarks,
                ref_fasta = ref_fasta,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_benchmark
        }
    }

    call Helpers.ConcatVcfs as MergeAnnotatedVcfs {
        input:
            vcfs = BenchmarkContig.annotated_vcf,
            vcfs_idx = BenchmarkContig.annotated_vcf_index,
            outfile_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }
    
    if (create_benchmarks) {
        call MergeSummaries {
            input:
                summary_files = select_all(BenchmarkContig.summary_file),
                prefix = prefix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_merge_summaries
        }

        call MergePlotTarballs {
            input:
                tarballs = select_all(BenchmarkContig.plot_tarball),
                prefix = prefix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_merge_tarballs
        }
    }


    output {
        File annotated_vcf = MergeAnnotatedVcfs.concat_vcf
        File annotated_vcf_index = MergeAnnotatedVcfs.concat_vcf_idx
        File? merged_summary = MergeSummaries.merged_file
        File? merged_plot_tarball = MergePlotTarballs.merged_tarball
    }
}

task BenchmarkContig {
    input {
        File vcf_eval
        File vcf_truth
        File sv_truth
        File closest_bed
        String contig
        String prefix
        Boolean create_benchmarks
        File ref_fasta
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String create_benchmarks_flag = if create_benchmarks then "--create_benchmarks" else ""

    command <<<
        set -euxo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/benchmark_annotations.py \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{sv_truth} \
            ~{closest_bed} \
            ~{contig} \
            ~{prefix} \
            --ref_fasta ~{ref_fasta} \
            ~{create_benchmarks_flag}
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.annotated.vcf.gz.tbi"
        File? summary_file = if create_benchmarks then "~{prefix}_benchmark_results/summary.txt" else ""
        File? plot_tarball = if create_benchmarks then "~{prefix}.benchmarks.tar.gz" else ""
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 32,
        disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB") + size(sv_truth, "GB")) * 3 + 50,
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

        mkdir final_results
        for tarball in ~{sep=' ' tarballs}; do
            tar -xzf $tarball -C final_results
        done

        tar -czf ~{prefix}.plots.tar.gz final_results/
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
