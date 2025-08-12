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

        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_subset_sv_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_annotate_benchmark
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_merge_benchmark_summaries
        RuntimeAttr? runtime_attr_merge_plot_tarballs

        RuntimeAttr? runtime_attr_truvari_filter_eval_vcf
        RuntimeAttr? runtime_attr_truvari_filter_truth_vcf
        RuntimeAttr? runtime_attr_truvari_run_truvari_09
        RuntimeAttr? runtime_attr_truvari_annotate_matched_09
        RuntimeAttr? runtime_attr_truvari_run_truvari_07
        RuntimeAttr? runtime_attr_truvari_annotate_matched_07
        RuntimeAttr? runtime_attr_truvari_run_truvari_05
        RuntimeAttr? runtime_attr_truvari_annotate_matched_05
        RuntimeAttr? runtime_attr_truvari_concat_matched

        RuntimeAttr? runtime_attr_bedtools_convert_to_symbolic
        RuntimeAttr? runtime_attr_bedtools_split_eval
        RuntimeAttr? runtime_attr_bedtools_split_truth
        RuntimeAttr? runtime_attr_bedtools_compare_del
        RuntimeAttr? runtime_attr_bedtools_calcu_del
        RuntimeAttr? runtime_attr_bedtools_compare_dup
        RuntimeAttr? runtime_attr_bedtools_calcu_dup
        RuntimeAttr? runtime_attr_bedtools_compare_ins
        RuntimeAttr? runtime_attr_bedtools_calcu_ins
        RuntimeAttr? runtime_attr_bedtools_compare_inv
        RuntimeAttr? runtime_attr_bedtools_calcu_inv
        RuntimeAttr? runtime_attr_bedtools_compare_bnd
        RuntimeAttr? runtime_attr_bedtools_calcu_bnd
        RuntimeAttr? runtime_attr_bedtools_merge_comparisons
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
                runtime_attr_override = runtime_attr_subset_eval
        }
        call Helpers.SubsetVcfToContig as SubsetTruth {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_index,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\"'",
                prefix = "~{prefix}.~{contig}.truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call Helpers.RenameVariantIds as RenameTruthIds {
            input:
                vcf = SubsetTruth.subset_vcf,
                vcf_index = SubsetTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.truth.renamed",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_rename_truth
        }

        call Helpers.SubsetVcfToContig as SubsetSVTruth {
            input:
                vcf = vcf_sv_truth,
                vcf_index = vcf_sv_truth_index,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'",
                prefix = "~{prefix}.~{contig}.sv_truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_sv_truth
        }

        call Helpers.RenameVariantIds as RenameSVTruthIds {
            input:
                vcf = SubsetSVTruth.subset_vcf,
                vcf_index = SubsetSVTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_rename_sv_truth
        }

        call ExactMatch {
            input:
                vcf_eval = SubsetEval.subset_vcf,
                vcf_truth = RenameTruthIds.renamed_vcf,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call Truvari.TruvariMatch as TruvariMatches {
            input:
                vcf_eval = ExactMatch.unmatched_vcf,
                vcf_eval_index = ExactMatch.unmatched_vcf_index,
                vcf_truth = RenameTruthIds.renamed_vcf,
                vcf_truth_index = RenameTruthIds.renamed_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                truvari_docker = truvari_docker,
                runtime_attr_filter_eval_vcf = runtime_attr_truvari_filter_eval_vcf,
                runtime_attr_filter_truth_vcf = runtime_attr_truvari_filter_truth_vcf,
                runtime_attr_run_truvari_09 = runtime_attr_truvari_run_truvari_09,
                runtime_attr_annotate_matched_09 = runtime_attr_truvari_annotate_matched_09,
                runtime_attr_run_truvari_07 = runtime_attr_truvari_run_truvari_07,
                runtime_attr_annotate_matched_07 = runtime_attr_truvari_annotate_matched_07,
                runtime_attr_run_truvari_05 = runtime_attr_truvari_run_truvari_05,
                runtime_attr_annotate_matched_05 = runtime_attr_truvari_annotate_matched_05,
                runtime_attr_concat_matched = runtime_attr_truvari_concat_matched
        }

        call Bedtools.BedtoolsClosestSV as BedtoolsClosest {
            input:
                vcf_eval = TruvariMatches.unmatched_vcf,
                vcf_eval_index = TruvariMatches.unmatched_vcf_index,
                vcf_sv_truth = RenameSVTruthIds.renamed_vcf,
                vcf_sv_truth_index = RenameSVTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}",
                sv_pipeline_docker = pipeline_docker,
                runtime_attr_convert_to_symbolic = runtime_attr_bedtools_convert_to_symbolic,
                runtime_attr_split_eval = runtime_attr_bedtools_split_eval,
                runtime_attr_split_truth = runtime_attr_bedtools_split_truth,
                runtime_attr_compare_del = runtime_attr_bedtools_compare_del,
                runtime_attr_calcu_del = runtime_attr_bedtools_calcu_del,
                runtime_attr_compare_dup = runtime_attr_bedtools_compare_dup,
                runtime_attr_calcu_dup = runtime_attr_bedtools_calcu_dup,
                runtime_attr_compare_ins = runtime_attr_bedtools_compare_ins,
                runtime_attr_calcu_ins = runtime_attr_bedtools_calcu_ins,
                runtime_attr_compare_inv = runtime_attr_bedtools_compare_inv,
                runtime_attr_calcu_inv = runtime_attr_bedtools_calcu_inv,
                runtime_attr_compare_bnd = runtime_attr_bedtools_compare_bnd,
                runtime_attr_calcu_bnd = runtime_attr_bedtools_calcu_bnd,
                runtime_attr_merge_comparisons = runtime_attr_bedtools_merge_comparisons
        }

        call AnnotateAndBenchmark {
            input:
                exact_matched_vcf = ExactMatch.matched_vcf,
                exact_matched_vcf_index = ExactMatch.matched_vcf_index,
                truvari_matched_vcf = TruvariMatches.matched_vcf,
                truvari_matched_vcf_index = TruvariMatches.matched_vcf_index,
                truvari_too_small_vcf = TruvariMatches.dropped_vcf,
                truvari_too_small_vcf_index = TruvariMatches.dropped_vcf_index,
                truvari_unmatched_vcf = TruvariMatches.unmatched_vcf,
                truvari_unmatched_vcf_index = TruvariMatches.unmatched_vcf_index,
                closest_bed = BedtoolsClosest.closest_bed,
                vcf_truth_snv = RenameTruthIds.renamed_vcf,
                vcf_truth_snv_index = RenameTruthIds.renamed_vcf_index,
                vcf_truth_sv = RenameSVTruthIds.renamed_vcf,
                vcf_truth_sv_index = RenameSVTruthIds.renamed_vcf_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
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

    call Helpers.MergeResults as MergeBenchmarkSummaries {
        input:
            tsvs = select_all(AnnotateAndBenchmark.benchmark_summary_tsv),
            merged_filename = "~{prefix}.benchmark_summary.tsv",
            hail_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call Helpers.MergeResults as MergeSummaryStats {
        input:
            tsvs = select_all(AnnotateAndBenchmark.summary_stats_tsv),
            merged_filename = "~{prefix}.summary_stats.tsv",
            hail_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call MergePlotTarballs {
        input:
            tarballs = select_all(AnnotateAndBenchmark.plot_tarball),
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_plot_tarballs
    }

    output {
        File annotated_vcf = MergeFinalVcfs.concat_vcf
        File annotated_vcf_index = MergeFinalVcfs.concat_vcf_idx
        File plots_tarball = MergePlotTarballs.merged_tarball
        File benchmark_summaries_tsv = MergeBenchmarkSummaries.merged_tsv
        File summary_stats_tsv = MergeSummaryStats.merged_tsv
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
        set -euo pipefail
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
        cpu_cores: 1, 
        mem_gb: 8, 
        disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) * 2 + 10,
        boot_disk_gb: 10, 
        preemptible_tries: 2, 
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

task AnnotateAndBenchmark {
    input {
        File exact_matched_vcf
        File exact_matched_vcf_index
        File truvari_matched_vcf
        File truvari_matched_vcf_index
        File truvari_too_small_vcf
        File truvari_too_small_vcf_index
        File truvari_unmatched_vcf
        File truvari_unmatched_vcf_index
        File closest_bed
        File vcf_truth_snv
        File vcf_truth_snv_index
        File vcf_truth_sv
        File vcf_truth_sv_index
        String contig
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/annotate_and_benchmark.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --exact_matched_vcf ~{exact_matched_vcf} \
            --truvari_matched_vcf ~{truvari_matched_vcf} \
            --truvari_too_small_vcf ~{truvari_too_small_vcf} \
            --truvari_unmatched_vcf ~{truvari_unmatched_vcf} \
            --closest_bed ~{closest_bed} \
            --vcf_truth_snv ~{vcf_truth_snv} \
            --vcf_truth_sv ~{vcf_truth_sv}
    >>>

    output {
        File final_vcf = "~{prefix}.final_annotated.vcf.gz"
        File final_vcf_index = "~{prefix}.final_annotated.vcf.gz.tbi"
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
        File plot_tarball = "~{prefix}.benchmarks.tar.gz"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 2, 
        mem_gb: 16, 
        disk_gb: 10 + ceil((size(exact_matched_vcf, "GB") + size(truvari_matched_vcf, "GB") + size(truvari_too_small_vcf, "GB") + size(truvari_unmatched_vcf, "GB") + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) * 1.5),
        boot_disk_gb: 15, 
        preemptible_tries: 2, 
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
        set -euo pipefail

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
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: ceil(size(tarballs, "GB")) * 2 + 10,
        boot_disk_gb: 10, 
        preemptible_tries: 2,
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
