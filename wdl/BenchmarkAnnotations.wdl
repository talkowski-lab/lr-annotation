version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers
import "BedtoolsClosestSV.wdl" as Bedtools
import "TruvariMatch.wdl" as Truvari

workflow BenchmarkAnnotations {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth
        File vcf_truth_index
        File vcf_sv_truth
        File vcf_sv_truth_index
        File ref_fasta
        File ref_fasta_fai
        
        File primary_contigs_list
        String pipeline_docker
        String truvari_docker
        String prefix
        Boolean create_benchmarks

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_truvari_match
        RuntimeAttr? runtime_attr_bedtools
        RuntimeAttr? runtime_attr_annotate_benchmark
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_merge_summaries
        RuntimeAttr? runtime_attr_merge_tarballs
    }

    Array[String] contigs = read_lines(primary_contigs_list)

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetEval {
            input:
                vcf = vcf_eval,
                vcf_index = vcf_eval_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}.eval",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset
        }
        call Helpers.SubsetVcfToContig as SubsetTruth {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}.truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset
        }
        call Helpers.SubsetVcfToContig as SubsetSVTruth {
            input:
                vcf = vcf_sv_truth,
                vcf_index = vcf_sv_truth_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv_truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call ExactMatch {
            input:
                vcf_eval = SubsetEval.subset_vcf,
                vcf_truth = SubsetTruth.subset_vcf,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call Truvari.TruvariMatch as TruvariMatches {
            input:
                vcf_eval_unmatched = ExactMatch.unmatched_vcf,
                vcf_eval_unmatched_index = ExactMatch.unmatched_vcf_index,
                vcf_truth = SubsetTruth.subset_vcf,
                vcf_truth_index = SubsetTruth.subset_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                truvari_docker = truvari_docker,
                runtime_attr_override = runtime_attr_truvari_match
        }

        call Bedtools.BedtoolsClosestSV as BedtoolsClosest {
            input:
                vcf_eval = TruvariMatches.unmatched_vcf,
                vcf_eval_index = TruvariMatches.unmatched_vcf_index,
                vcf_sv_truth = SubsetSVTruth.subset_vcf,
                vcf_sv_truth_index = SubsetSVTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}",
                sv_pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_bedtools
        }

        call AnnotateAndBenchmark {
            input:
                vcf_unmatched_from_truvari = TruvariMatches.unmatched_vcf,
                closest_bed = BedtoolsClosest.closest_bed,
                exact_matched_vcf = ExactMatch.matched_vcf,
                truvari_matched_vcf = TruvariMatches.matched_vcf,
                vcf_truth_snv = SubsetTruth.subset_vcf,
                vcf_truth_sv = SubsetSVTruth.subset_vcf,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                create_benchmarks = create_benchmarks,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_annotate_benchmark
        }
    }

    call Helpers.ConcatVcfs as MergeFinalVcfs {
        input:
            vcfs = select_all(AnnotateAndBenchmark.final_vcf),
            vcfs_idx = select_all(AnnotateAndBenchmark.final_vcf_index),
            outfile_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }
    
    if (create_benchmarks) {
        call MergeSummaries {
            input:
                summary_files = select_all(AnnotateAndBenchmark.summary_file),
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_summaries
    }

    call MergePlotTarballs {
        input:
                tarballs = select_all(AnnotateAndBenchmark.plot_tarball),
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_tarballs
        }
    }

    output {
        File annotated_vcf = MergeFinalVcfs.concat_vcf
        File annotated_vcf_index = MergeFinalVcfs.concat_vcf_idx
        File? merged_summary = MergeSummaries.merged_file
        File? merged_plot_tarball = MergePlotTarballs.merged_tarball
    }
}

task ExactMatch {
    input {
        File vcf_eval
        File vcf_truth
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/exact_match.py \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{prefix}
    >>>

    output {
        File matched_vcf = "~{prefix}.exact_matched.vcf.gz"
        File matched_vcf_index = "~{prefix}.exact_matched.vcf.gz.tbi"
        File unmatched_vcf = "~{prefix}.unmatched.vcf.gz"
        File unmatched_vcf_index = "~{prefix}.unmatched.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 8, disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) * 2 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
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

task AnnotateAndBenchmark {
    input {
        File vcf_unmatched_from_truvari
        File closest_bed
        File exact_matched_vcf
        File truvari_matched_vcf
        File vcf_truth_snv
        File vcf_truth_sv
        String contig
        String prefix
        Boolean create_benchmarks
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String create_benchmarks_flag = if create_benchmarks then "--create_benchmarks" else ""
    String vcf_truth_snv_arg = if create_benchmarks then "--vcf_truth_snv ~{vcf_truth_snv}" else ""
    String vcf_truth_sv_arg = if create_benchmarks then "--vcf_truth_sv ~{vcf_truth_sv}" else ""
    String exact_matched_arg = if create_benchmarks then "--exact_matched_vcf ~{exact_matched_vcf}" else ""
    String truvari_matched_arg = if create_benchmarks then "--truvari_matched_vcf ~{truvari_matched_vcf}" else ""
    String contig_arg = if create_benchmarks then "--contig ~{contig}" else ""


    command <<<
        set -euxo pipefail

        # Run the main python script to annotate bedtools matches and optionally create benchmarks
        python3 /opt/gnomad-lr/scripts/benchmark/annotate_and_benchmark.py \
            ~{vcf_unmatched_from_truvari} \
            ~{closest_bed} \
            ~{prefix} \
            ~{create_benchmarks_flag} \
            ~{vcf_truth_snv_arg} \
            ~{vcf_truth_sv_arg} \
            ~{exact_matched_arg} \
            ~{truvari_matched_arg} \
            ~{contig_arg}

        # Combine all parts into the final VCF for the contig
        bcftools concat -a -f -Oz -o ~{prefix}.final_annotated.vcf.gz \
            ~{exact_matched_vcf} \
            ~{truvari_matched_vcf} \
            ~{prefix}.bedtools_matched.vcf.gz \
            ~{prefix}.final_unmatched.vcf.gz
        
        tabix -p vcf -f ~{prefix}.final_annotated.vcf.gz
    >>>

    output {
        File final_vcf = "~{prefix}.final_annotated.vcf.gz"
        File final_vcf_index = "~{prefix}.final_annotated.vcf.gz.tbi"
        File? summary_file = if create_benchmarks then "~{prefix}_benchmark_results/summary.txt" else ""
        File? plot_tarball = if create_benchmarks then "~{prefix}.benchmarks.tar.gz" else ""
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 2, mem_gb: 16, disk_gb: ceil(size(vcf_unmatched_from_truvari, "GB") * 2 + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB") + size(exact_matched_vcf, "GB") + size(truvari_matched_vcf, "GB")) + 30,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
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
        cpu_cores: 1, mem_gb: 2, disk_gb: ceil(size(summary_files, "GB")) * 2 + 5,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
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

        mkdir -p final_results/AF_plots
        mkdir -p final_results/VEP_plots

        for tarball in ~{sep=' ' tarballs}; do
            tar -xvf $tarball --strip-components=1 -C final_results
        done

        tar -czf ~{prefix}.plots.tar.gz final_results/
    >>>

    output {
        File merged_tarball = "~{prefix}.plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 4, disk_gb: ceil(size(tarballs, "GB")) * 2 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
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
