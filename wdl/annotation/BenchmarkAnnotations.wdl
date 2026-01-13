version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers
import "../benchmarking/BedtoolsClosestSV.wdl"
import "../benchmarking/TruvariMatch.wdl"

workflow BenchmarkAnnotations {
    input {
        File vcf
        File vcf_idx
        File vcf_truth
        File vcf_truth_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx

        String prefix
        Int variants_per_shard
        Int truvari_match_min_length = 10
        String? skip_vep_categories = "hgvsc,cdna_position,distance,hgvsp,domains,ensp"

        File ref_fa
        File ref_fai
        Array[String] contigs
        
        String benchmark_annotations_docker
        String utils_docker

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
        RuntimeAttr? runtime_attr_collect_matched_ids
        RuntimeAttr? runtime_attr_extract_eval_vep_header
        RuntimeAttr? runtime_attr_extract_truth_vep_header
        RuntimeAttr? runtime_attr_shard_matched_eval
        RuntimeAttr? runtime_attr_compute_shard_benchmarks
        RuntimeAttr? runtime_attr_merge_shard_benchmarks
        RuntimeAttr? runtime_attr_compute_summary_for_contig
        RuntimeAttr? runtime_attr_merge_benchmark_summaries
        RuntimeAttr? runtime_attr_merge_plot_tarballs
    }

    scatter (contig in contigs) {
        call Helpers.StripGenotypes as StripEvalGenotypes {
            input:
                vcf = vcf,
                vcf_index = vcf_idx,
                prefix = "~{prefix}.~{contig}.eval.stripped",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_eval
        }

        call Helpers.StripGenotypes as StripTruthGenotypes {
            input:
                vcf = vcf_truth,
                vcf_index = vcf_truth_idx,
                prefix = "~{prefix}.~{contig}.truth.stripped",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call Helpers.StripGenotypes as StripSVTruthGenotypes {
            input:
                vcf = vcf_sv_truth,
                vcf_index = vcf_sv_truth_idx,
                prefix = "~{prefix}.~{contig}.sv_truth.stripped",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_sv_truth
        }

        call Helpers.SubsetVcfToContig as SubsetEval {
            input:
                vcf = StripEvalGenotypes.stripped_vcf,
                vcf_index = StripEvalGenotypes.stripped_vcf_index,
                contig = contig,
                strip_genotypes = false,
                prefix = "~{prefix}.~{contig}.eval",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_eval
        }

        call Helpers.SubsetVcfToContig as SubsetTruth {
            input:
                vcf = StripTruthGenotypes.stripped_vcf,
                vcf_index = StripTruthGenotypes.stripped_vcf_index,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\"'",
                strip_genotypes = false,
                prefix = "~{prefix}.~{contig}.truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call Helpers.SubsetVcfToContig as SubsetSVTruth {
            input:
                vcf = StripSVTruthGenotypes.stripped_vcf,
                vcf_index = StripSVTruthGenotypes.stripped_vcf_index,
                contig = contig,
                args_string = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'",
                strip_genotypes = false,
                prefix = "~{prefix}.~{contig}.sv_truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_sv_truth
        }

        call Helpers.RenameVariantIds as RenameTruthIds {
            input:
                vcf = SubsetTruth.subset_vcf,
                vcf_index = SubsetTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.truth.renamed",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_rename_truth
        }

        call Helpers.RenameVariantIds as RenameSVTruthIds {
            input:
                vcf = SubsetSVTruth.subset_vcf,
                vcf_index = SubsetSVTruth.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_rename_sv_truth
        }

        call ExtractVepHeader as ExtractTruthVepHeader {
            input:
                vcf = RenameTruthIds.renamed_vcf,
                vcf_index = RenameTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}.truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_extract_truth_vep_header
        }

        call ExtractVepHeader as ExtractEvalVepHeader {
            input:
                vcf = SubsetEval.subset_vcf,
                vcf_index = SubsetEval.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.eval",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_extract_eval_vep_header
        }

        call ExactMatch {
            input:
                vcf_eval = SubsetEval.subset_vcf,
                vcf_truth = RenameTruthIds.renamed_vcf,
                prefix = "~{prefix}.~{contig}",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call TruvariMatch.TruvariMatch {
            input:
                vcf_eval = ExactMatch.unmatched_vcf,
                vcf_eval_index = ExactMatch.unmatched_vcf_index,
                vcf_truth = RenameTruthIds.renamed_vcf,
                vcf_truth_index = RenameTruthIds.renamed_vcf_index,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}",
                min_sv_length = truvari_match_min_length,
                utils_docker = utils_docker,
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

        call BedtoolsClosestSV.BedtoolsClosestSV {
            input:
                vcf_eval = TruvariMatch.unmatched_vcf,
                vcf_eval_index = TruvariMatch.unmatched_vcf_index,
                vcf_sv_truth = RenameSVTruthIds.renamed_vcf,
                vcf_sv_truth_index = RenameSVTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}",
                bedtools_closest_docker = benchmark_annotations_docker,
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

        call Helpers.ConcatTsvs as BuildAnnotationTsv {
            input:
                tsvs = [ExactMatch.annotation_tsv, TruvariMatch.annotation_tsv, BedtoolsClosestSV.annotation_tsv],
                prefix = "~{prefix}.~{contig}.annotations",
                docker = benchmark_annotations_docker
        }

        call CollectMatchedIDsAndINFO {
            input:
                annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                vcf_eval = SubsetEval.subset_vcf,
                vcf_eval_index = SubsetEval.subset_vcf_index,
                vcf_truth_snv = RenameTruthIds.renamed_vcf,
                vcf_truth_snv_index = RenameTruthIds.renamed_vcf_index,
                vcf_truth_sv = RenameSVTruthIds.renamed_vcf,
                vcf_truth_sv_index = RenameSVTruthIds.renamed_vcf_index,
                prefix = "~{prefix}.~{contig}",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_collect_matched_ids
        }

        call ShardedMatchedVariants {
            input:
                matched_with_info_tsv = CollectMatchedIDsAndINFO.matched_with_info_tsv,
                variants_per_shard = variants_per_shard,
                prefix = "~{prefix}.~{contig}",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_shard_matched_eval
        }

        scatter (shard_idx in range(length(ShardedMatchedVariants.shard_tsvs))) {
            call ComputeShardBenchmarks {
                input:
                    matched_shard_tsv = ShardedMatchedVariants.shard_tsvs[shard_idx],
                    eval_vep_header = ExtractEvalVepHeader.vep_header_txt,
                    truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                    contig = contig,
                    shard_label = "~{shard_idx}",
                    prefix = "~{prefix}.~{contig}",
                    skip_vep_categories = skip_vep_categories,
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_compute_shard_benchmarks
            }
        }

        call MergeShardBenchmarks {
            input:
                af_pair_tsvs = select_all(ComputeShardBenchmarks.af_pairs_tsv),
                vep_pair_tsvs = select_all(ComputeShardBenchmarks.vep_pairs_tsv),
                truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                skip_vep_categories = skip_vep_categories,
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_merge_shard_benchmarks
        }

        call ComputeSummaryForContig {
            input:
                final_vcf = SubsetEval.subset_vcf,
                final_vcf_index = SubsetEval.subset_vcf_index,
                annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                matched_with_info_tsv = CollectMatchedIDsAndINFO.matched_with_info_tsv,
                truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_compute_summary_for_contig
        }
    }

    call Helpers.ConcatTsvs as MergeAnnotationTsvs {
        input:
            tsvs = select_all(BuildAnnotationTsv.concatenated_tsv),
            prefix = "~{prefix}.annotations",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call Helpers.ConcatTsvs as MergeBenchmarkSummaries {
        input:
            tsvs = select_all(ComputeSummaryForContig.benchmark_summary_tsv),
            prefix = "~{prefix}.benchmark_summary",
            preserve_header = true,
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call Helpers.ConcatTsvs as MergeSummaryStats {
        input:
            tsvs = select_all(ComputeSummaryForContig.summary_stats_tsv),
            prefix = "~{prefix}.summary_stats",
            preserve_header = true,
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_merge_benchmark_summaries
    }

    call MergePlotTarballs {
        input:
            tarballs = select_all(MergeShardBenchmarks.plot_tarball),
            prefix = prefix,
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_merge_plot_tarballs
    }

    output {
        File annotations_tsv_benchmark = MergeAnnotationTsvs.concatenated_tsv
        File plots_tarball = MergePlotTarballs.merged_tarball
        File benchmark_summaries_tsv = MergeBenchmarkSummaries.concatenated_tsv
        File summary_stats_tsv = MergeSummaryStats.concatenated_tsv
    }
}

task ExactMatch {
    input {
        File vcf_eval
        File vcf_truth
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/exact_match.py \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{prefix}
        
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/gnomAD_V4_match\t%INFO/gnomAD_V4_match_ID\n' \
            ~{prefix}.exact_matched.vcf.gz > ~{prefix}.exact_matched.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.exact_matched.tsv"
        File unmatched_vcf = "~{prefix}.unmatched.vcf.gz"
        File unmatched_vcf_index = "~{prefix}.unmatched.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) * 1.5 + 5,
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

task CollectMatchedIDsAndINFO {
    input {
        File annotation_tsv
        File vcf_eval
        File vcf_eval_index
        File vcf_truth_snv
        File vcf_truth_snv_index
        File vcf_truth_sv
        File vcf_truth_sv_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -F'\t' '{print $5}' ~{annotation_tsv} | sort -u > eval_ids.list
        awk -F'\t' '{print $7}' ~{annotation_tsv} | sort -u > truth_ids.list

        bcftools view -i 'ID=@eval_ids.list' ~{vcf_eval} \
            | bcftools query -f '%ID\t%INFO\n' \
            | sort -k1,1 > eval_info.tsv

        bcftools view -i 'ID=@truth_ids.list' ~{vcf_truth_snv} \
            | bcftools query -f '%ID\t%INFO\n' > truth_snv_info.tsv
        bcftools view -i 'ID=@truth_ids.list' ~{vcf_truth_sv} \
            | bcftools query -f '%ID\t%INFO\n' > truth_sv_info.tsv
        cat truth_snv_info.tsv truth_sv_info.tsv | sort -k1,1 > truth_info.tsv

        awk -F'\t' 'BEGIN{OFS="\t"} {print $5"\t"$7}' ~{annotation_tsv} \
            | sort -k1,1 \
            | join -t $'\t' -1 1 -2 1 - eval_info.tsv \
            | sort -k2,2 \
            | join -t $'\t' -1 2 -2 1 - truth_info.tsv \
            | awk 'BEGIN{OFS="\t"} {print $2,$1,$3,$4}' \
            | bgzip -c > ~{prefix}.matched_with_info.tsv.gz
        
        tabix -s 1 -b 1 -e 1 ~{prefix}.matched_with_info.tsv.gz
    >>>

    output {
        File matched_with_info_tsv = "~{prefix}.matched_with_info.tsv.gz"
        File matched_with_info_tsv_tbi = "~{prefix}.matched_with_info.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) + 10,
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

task ExtractVepHeader {
    input {
        File vcf
        File vcf_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools view \
            -h \
            ~{vcf} \
        | awk 'BEGIN{IGNORECASE=1} /^##INFO=<ID=(vep|csq),/ {print; exit}' > ~{prefix}_vep_header.txt
    >>>

    output {
        File vep_header_txt = "~{prefix}_vep_header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 5,
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

task ShardedMatchedVariants {
    input {
        File matched_with_info_tsv
        Int variants_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        mkdir -p shards
        
        if [ ! -s ~{matched_with_info_tsv} ] || [ $(zcat ~{matched_with_info_tsv} | wc -l) -eq 0 ]; then
            echo "No matched variants found, creating empty shard"
            touch shards/matched.000000.tsv
            bgzip -f shards/matched.000000.tsv
        else
            zcat ~{matched_with_info_tsv} | awk 'BEGIN{c=0;f=0} {print > sprintf("shards/matched.%06d.tsv", int(c/~{variants_per_shard})) ; c++} END{ }'
            if ls shards/matched.*.tsv 1> /dev/null 2>&1; then
                ls shards/matched.*.tsv | while read f; do bgzip -f "$f"; done
            fi
        fi
    >>>

    output {
        Array[File] shard_tsvs = glob("shards/*.tsv.gz")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(matched_with_info_tsv, "GB")) * 2 + 5,
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

task ComputeShardBenchmarks {
    input {
        File matched_shard_tsv
        File eval_vep_header
        File truth_vep_header
        String contig
        String shard_label
        String prefix
        String? skip_vep_categories
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ ! -s ~{matched_shard_tsv} ] || [ $(zcat ~{matched_shard_tsv} | wc -l) -eq 0 ]; then
            echo "Empty shard, creating empty output files"
            echo -e "af_key\teval_af\ttruth_af" | gzip > ~{prefix}.shard_~{shard_label}.af_pairs.tsv.gz
            echo -e "category\teval\ttruth\tcount" | gzip > ~{prefix}.shard_~{shard_label}.vep_pairs.tsv.gz
        else
            python3 /opt/gnomad-lr/scripts/benchmark/compute_benchmarks_shard.py \
                --prefix ~{prefix} \
                --contig ~{contig} \
                --matched_shard_tsv ~{matched_shard_tsv} \
                --eval_vep_header ~{eval_vep_header} \
                --truth_vep_header ~{truth_vep_header} \
                --shard_label ~{shard_label} \
                ~{"--skip_vep_categories " + skip_vep_categories}
        fi
    >>>

    output {
        File af_pairs_tsv = "~{prefix}.shard_~{shard_label}.af_pairs.tsv.gz"
        File vep_pairs_tsv = "~{prefix}.shard_~{shard_label}.vep_pairs.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 5,
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

task MergeShardBenchmarks {
    input {
        Array[File] af_pair_tsvs
        Array[File] vep_pair_tsvs
        File truth_vep_header
        String contig
        String prefix
        String? skip_vep_categories
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/merge_benchmarks_from_pairs.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --af_pair_tsvs ~{sep=',' af_pair_tsvs} \
            --vep_pair_tsvs ~{sep=',' vep_pair_tsvs} \
            --truth_vep_header ~{truth_vep_header} \
            ~{"--skip_vep_categories " + skip_vep_categories}
    >>>

    output {
        File plot_tarball = "~{prefix}.benchmarks.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
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

task ComputeSummaryForContig {
    input {
        File final_vcf
        File final_vcf_index
        File annotation_tsv
        File matched_with_info_tsv
        File truth_vep_header
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_summary_for_contig.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --final_vcf ~{final_vcf} \
            --annotation_tsv ~{annotation_tsv} \
            --matched_with_info_tsv ~{matched_with_info_tsv} \
            --truth_vep_header ~{truth_vep_header}
    >>>

    output {
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(final_vcf, "GB")) + 5,
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

task MergePlotTarballs {
    input {
        Array[File] tarballs
        String prefix
        String docker
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
        disk_gb: ceil(size(tarballs, "GB") * 1.5) + 5,
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
