version 1.0

import "general/Structs.wdl"
import "general/Helpers.wdl" as Helpers
import "benchmarking/BedtoolsClosestSV.wdl" as Bedtools
import "benchmarking/TruvariMatch.wdl" as Truvari
import "benchmarking/ShardedBenchmarks.wdl" as Sharded

workflow BenchmarkAnnotations {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx

        File ref_fasta
        File ref_fasta_fai
        File primary_contigs_list
        
        String pipeline_docker
        String truvari_docker
        String prefix

        Int variants_per_shard
        String? skip_vep_categories = "hgvsc,cdna_position,distance,hgvsp,domains,ensp"

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

        RuntimeAttr? runtime_attr_bedtools_annotate_unmatched
        RuntimeAttr? runtime_attr_build_final_vcf
        RuntimeAttr? runtime_attr_collect_matched_ids
        RuntimeAttr? runtime_attr_extract_truth_vep_header
        RuntimeAttr? runtime_attr_extract_truth_info
        RuntimeAttr? runtime_attr_shard_matched_eval
        RuntimeAttr? runtime_attr_compute_shard_benchmarks
        RuntimeAttr? runtime_attr_merge_shard_benchmarks
        RuntimeAttr? runtime_attr_compute_summary_for_contig
    }

    Array[String] contigs = read_lines(primary_contigs_list)

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetEval {
            input:
                vcf = vcf_eval,
                vcf_index = vcf_eval_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.eval",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_eval
        }

        call Helpers.SubsetVcfToContig as SubsetTruth {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_idx,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\"'",
                prefix = "~{prefix}.~{contig}.truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call Helpers.SubsetVcfToContig as SubsetSVTruth {
            input:
                vcf = vcf_sv_truth,
                vcf_index = vcf_sv_truth_idx,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'",
                prefix = "~{prefix}.~{contig}.sv_truth",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_subset_sv_truth
        }

        call Helpers.RenameVariantIds as RenameTruthIds {
            input:
                vcf = SubsetTruth.subset_vcf,
                vcf_index = SubsetTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.truth.renamed",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_rename_truth
        }

        call Helpers.RenameVariantIds as RenameSVTruthIds {
            input:
                vcf = SubsetSVTruth.subset_vcf,
                vcf_index = SubsetSVTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_rename_sv_truth
        }

        call ExtractTruthVepHeader {
            input:
                vcf_truth_snv = RenameTruthIds.renamed_vcf,
                vcf_truth_snv_index = RenameTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_extract_truth_vep_header
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

        call AnnotateBedtoolsMatches {
            input:
                truvari_unmatched_vcf = TruvariMatches.unmatched_vcf,
                truvari_unmatched_vcf_index = TruvariMatches.unmatched_vcf_index,
                closest_bed = BedtoolsClosest.closest_bed,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_bedtools_annotate_unmatched
        }

        call Helpers.ConcatVcfs as BuildFinalContigVcf {
            input:
                vcfs = [ExactMatch.matched_vcf, TruvariMatches.matched_vcf, TruvariMatches.dropped_vcf, AnnotateBedtoolsMatches.bedtools_matched_vcf, AnnotateBedtoolsMatches.final_unmatched_vcf],
                vcfs_idx = [ExactMatch.matched_vcf_index, TruvariMatches.matched_vcf_index, TruvariMatches.dropped_vcf_index, AnnotateBedtoolsMatches.bedtools_matched_vcf_index, AnnotateBedtoolsMatches.final_unmatched_vcf_index],
                outfile_prefix = "~{prefix}.~{contig}",
                docker_image = pipeline_docker,
                runtime_attr_override = runtime_attr_build_final_vcf
        }

        call CollectMatchedIDs {
            input:
                final_vcf = BuildFinalContigVcf.concat_vcf,
                final_vcf_index = BuildFinalContigVcf.concat_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_collect_matched_ids
        }

        call ExtractTruthInfoForMatched {
            input:
                matched_ids_tsv = CollectMatchedIDs.matched_ids_tsv,
                vcf_truth_snv = RenameTruthIds.renamed_vcf,
                vcf_truth_snv_index = RenameTruthIds.renamed_vcf_index,
                vcf_truth_sv = RenameSVTruthIds.renamed_vcf,
                vcf_truth_sv_index = RenameSVTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_extract_truth_info
        }

        call Sharded.ShardAndComputeBenchmarks as RunShardedBenchmarks {
            input:
                final_vcf = BuildFinalContigVcf.concat_vcf,
                final_vcf_index = BuildFinalContigVcf.concat_vcf_idx,
                matched_ids_tsv = CollectMatchedIDs.matched_ids_tsv,
                truth_tsv_snv = ExtractTruthInfoForMatched.truth_info_snv_tsv,
                truth_tsv_sv = ExtractTruthInfoForMatched.truth_info_sv_tsv,
                truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                variants_per_shard = variants_per_shard,
                skip_vep_categories = skip_vep_categories,
                runtime_attr_shard_matched_eval = runtime_attr_shard_matched_eval,
                runtime_attr_compute_shard_benchmarks = runtime_attr_compute_shard_benchmarks,
                runtime_attr_merge_shard_benchmarks = runtime_attr_merge_shard_benchmarks
        }

        call ComputeSummaryForContig {
            input:
                final_vcf = BuildFinalContigVcf.concat_vcf,
                final_vcf_index = BuildFinalContigVcf.concat_vcf_idx,
                truth_tsv_snv = ExtractTruthInfoForMatched.truth_info_snv_tsv,
                truth_tsv_sv = ExtractTruthInfoForMatched.truth_info_sv_tsv,
                truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_compute_summary_for_contig
        }
    }

    call Helpers.ConcatVcfs as MergeFinalVcfs {
        input:
            vcfs = select_all(BuildFinalContigVcf.concat_vcf),
            vcfs_idx = select_all(BuildFinalContigVcf.concat_vcf_idx),
            outfile_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    call Helpers.MergeResults as MergeBenchmarkSummaries {
        input:
            tsvs = select_all(ComputeSummaryForContig.benchmark_summary_tsv),
            merged_filename = "~{prefix}.benchmark_summary.tsv",
            hail_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call Helpers.MergeResults as MergeSummaryStats {
        input:
            tsvs = select_all(ComputeSummaryForContig.summary_stats_tsv),
            merged_filename = "~{prefix}.summary_stats.tsv",
            hail_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call MergePlotTarballs {
        input:
            tarballs = select_all(RunShardedBenchmarks.plot_tarball),
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

task AnnotateBedtoolsMatches {
    input {
        File truvari_unmatched_vcf
        File truvari_unmatched_vcf_index
        File closest_bed
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=gnomAD_V4_match,Number=1,Type=String,Description="Matching status against gnomAD v4.">' > header.hdr
        echo '##INFO=<ID=gnomAD_V4_match_ID,Number=1,Type=String,Description="Matching variant ID from gnomAD v4.">' >> header.hdr

        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' ~{truvari_unmatched_vcf} | sort -k1,1 > eval.coords.tsv
        awk '($1 != "query_svid"){print $1"\t"$2}' ~{closest_bed} | sort -k1,1 > closest.map.tsv
        join -t $'\t' -1 1 -2 1 closest.map.tsv eval.coords.tsv | awk -F'\t' '($2 != "."){print $3"\t"$4"\t"$5"\t"$6"\tBEDTOOLS_CLOSEST\t"$2}' | bgzip -c > annots.tab.gz
        tabix -s 1 -b 2 -e 2 annots.tab.gz

        bcftools annotate -a annots.tab.gz -h header.hdr -c CHROM,POS,REF,ALT,gnomAD_V4_match,gnomAD_V4_match_ID -Oz -o ~{prefix}.unmatched.annotated.vcf.gz ~{truvari_unmatched_vcf}
        tabix -p vcf -f ~{prefix}.unmatched.annotated.vcf.gz

        bcftools view -i 'INFO/gnomAD_V4_match!="."' ~{prefix}.unmatched.annotated.vcf.gz -Oz -o ~{prefix}.bedtools_matched.vcf.gz
        tabix -p vcf -f ~{prefix}.bedtools_matched.vcf.gz

        bcftools view -e 'INFO/gnomAD_V4_match!="."' ~{prefix}.unmatched.annotated.vcf.gz -Oz -o ~{prefix}.final_unmatched.vcf.gz
        tabix -p vcf -f ~{prefix}.final_unmatched.vcf.gz
    >>>

    output {
        File bedtools_matched_vcf = "~{prefix}.bedtools_matched.vcf.gz"
        File bedtools_matched_vcf_index = "~{prefix}.bedtools_matched.vcf.gz.tbi"
        File final_unmatched_vcf = "~{prefix}.final_unmatched.vcf.gz"
        File final_unmatched_vcf_index = "~{prefix}.final_unmatched.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10 + ceil(size(truvari_unmatched_vcf, "GB")) * 3,
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


task CollectMatchedIDs {
    input {
        File final_vcf
        File final_vcf_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        bcftools query -i 'INFO/gnomAD_V4_match!="."' -f '%ID\t%INFO/gnomAD_V4_match_ID\n' ~{final_vcf} | bgzip -c > ~{prefix}.matched_ids.tsv.gz
        tabix -s 1 -b 1 -e 1 ~{prefix}.matched_ids.tsv.gz
    >>>

    output {
        File matched_ids_tsv = "~{prefix}.matched_ids.tsv.gz"
        File matched_ids_tsv_tbi = "~{prefix}.matched_ids.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(final_vcf, "GB")) * 2 + 5,
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


task ExtractTruthVepHeader {
    input {
        File vcf_truth_snv
        File vcf_truth_snv_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools view -h ~{vcf_truth_snv} | awk 'BEGIN{IGNORECASE=1} /^##INFO=<ID=(vep|csq),/ {print; exit}' > ~{prefix}.truth_vep_header.txt
    >>>

    output {
        File vep_header_txt = "~{prefix}.truth_vep_header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf_truth_snv, "GB") * 1.5 + 2),
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


task ExtractTruthInfoForMatched {
    input {
        File matched_ids_tsv
        File vcf_truth_snv
        File vcf_truth_snv_index
        File vcf_truth_sv
        File vcf_truth_sv_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        zcat ~{matched_ids_tsv} | cut -f2 | sort -u > match_ids.list
        bcftools view -i 'ID=@match_ids.list' ~{vcf_truth_snv} | bcftools query -f '%ID\t%INFO\n' | bgzip -c > ~{prefix}.truth_info_snv.tsv.gz
        tabix -s 1 -b 1 -e 1 ~{prefix}.truth_info_snv.tsv.gz
        bcftools view -i 'ID=@match_ids.list' ~{vcf_truth_sv} | bcftools query -f '%ID\t%INFO\n' | bgzip -c > ~{prefix}.truth_info_sv.tsv.gz
        tabix -s 1 -b 1 -e 1 ~{prefix}.truth_info_sv.tsv.gz
    >>>

    output {
        File truth_info_snv_tsv = "~{prefix}.truth_info_snv.tsv.gz"
        File truth_info_snv_tsv_tbi = "~{prefix}.truth_info_snv.tsv.gz.tbi"
        File truth_info_sv_tsv = "~{prefix}.truth_info_sv.tsv.gz"
        File truth_info_sv_tsv_tbi = "~{prefix}.truth_info_sv.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) * 2 + 10,
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

task ComputeSummaryForContig {
    input {
        File final_vcf
        File final_vcf_index
        File truth_tsv_snv
        File truth_tsv_sv
        File truth_vep_header
        String contig
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/compute_summary_for_contig.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --final_vcf ~{final_vcf} \
            --truth_tsv_snv ~{truth_tsv_snv} \
            --truth_tsv_sv ~{truth_tsv_sv} \
            --truth_vep_header ~{truth_vep_header}
    >>>

    output {
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(final_vcf, "GB")) * 2 + 5,
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
