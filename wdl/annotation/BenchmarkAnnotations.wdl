version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers
import "../benchmarking/BedtoolsClosestSV.wdl"
import "../benchmarking/TruvariMatch.wdl"

workflow BenchmarkAnnotations {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx
        Array[String] contigs
        String prefix
        
        Int variants_per_shard
        Int min_svlen_eval
        Int min_svlen_truth
        String? skip_vep_categories = "hgvsc,cdna_position,distance,hgvsp,domains,ensp"
        
        String? args_string_vcf
        String? args_string_vcf_truth
        String? args_string_vcf_sv_truth
        String? rename_id_string_vcf
        String? rename_id_string_vcf_truth
        String? rename_id_string_vcf_sv_truth
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_vcf_truth
        Boolean? rename_id_strip_chr_vcf_sv_truth

        File ref_fa
        File ref_fai
        
        String benchmark_annotations_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_subset_sv_truth
        RuntimeAttr? runtime_attr_rename_eval
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_annotate_benchmark
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_merge_benchmark_summaries
        RuntimeAttr? runtime_attr_merge_plot_tarballs

        RuntimeAttr? runtime_attr_truvari_subset_eval
        RuntimeAttr? runtime_attr_truvari_subset_truth
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
        call Helpers.SubsetVcfToContig as SubsetEval {
            input:
                vcf = vcf_eval,
                vcf_idx = vcf_eval_idx,
                contig = contig,
                args_string = args_string_vcf,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.eval",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_eval
        }

        call Helpers.SubsetVcfToContig as SubsetTruth {
            input:
                vcf = vcf_truth,
                vcf_idx = vcf_truth_idx,
                contig = contig,
                args_string = args_string_vcf_truth,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call Helpers.SubsetVcfToContig as SubsetSVTruth {
            input:
                vcf = vcf_sv_truth,
                vcf_idx = vcf_sv_truth_idx,
                contig = contig,
                args_string = args_string_vcf_sv_truth,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.sv_truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_subset_sv_truth
        }

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIds {
                input:
                    vcf = SubsetEval.subset_vcf,
                    vcf_idx = SubsetEval.subset_vcf_idx,
                    id_format = select_first([rename_id_string_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                    prefix = "~{prefix}.~{contig}.eval.renamed",
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_rename_eval
            }
        }

        if (defined(rename_id_string_vcf_truth)) {
            call Helpers.RenameVariantIds as RenameTruthIds {
                input:
                    vcf = SubsetTruth.subset_vcf,
                    vcf_idx = SubsetTruth.subset_vcf_idx,
                    id_format = select_first([rename_id_string_vcf_truth]),
                    strip_chr = select_first([rename_id_strip_chr_vcf_truth, false]),
                    prefix = "~{prefix}.~{contig}.truth.renamed",
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_rename_truth
            }
        }

        if (defined(rename_id_string_vcf_sv_truth)) {
            call Helpers.RenameVariantIds as RenameSVTruthIds {
                input:
                    vcf = SubsetSVTruth.subset_vcf,
                    vcf_idx = SubsetSVTruth.subset_vcf_idx,
                    id_format = select_first([rename_id_string_vcf_sv_truth]),
                    strip_chr = select_first([rename_id_strip_chr_vcf_sv_truth, false]),
                    prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_rename_sv_truth
            }
        }

        File eval_vcf_final = select_first([RenameEvalIds.renamed_vcf, SubsetEval.subset_vcf])
        File eval_vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, SubsetEval.subset_vcf_idx])
        File truth_vcf_final = select_first([RenameTruthIds.renamed_vcf, SubsetTruth.subset_vcf])
        File truth_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, SubsetTruth.subset_vcf_idx])
        File sv_truth_vcf_final = select_first([RenameSVTruthIds.renamed_vcf, SubsetSVTruth.subset_vcf])
        File sv_truth_vcf_final_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, SubsetSVTruth.subset_vcf_idx])

        call ExtractVepHeader as ExtractTruthVepHeader {
            input:
                vcf = truth_vcf_final,
                vcf_idx = truth_vcf_final_idx,
                prefix = "~{prefix}.~{contig}.truth",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_extract_truth_vep_header
        }

        call ExtractVepHeader as ExtractEvalVepHeader {
            input:
                vcf = eval_vcf_final,
                vcf_idx = eval_vcf_final_idx,
                prefix = "~{prefix}.~{contig}.eval",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_extract_eval_vep_header
        }

        call ExactMatch {
            input:
                vcf_eval = eval_vcf_final,
                vcf_eval_idx = eval_vcf_final_idx,
                vcf_truth = truth_vcf_final,
                vcf_truth_idx = truth_vcf_final_idx,
                prefix = "~{prefix}.~{contig}",
                docker = benchmark_annotations_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call TruvariMatch.TruvariMatch {
            input:
                vcf_eval = ExactMatch.unmatched_vcf,
                vcf_eval_idx = ExactMatch.unmatched_vcf_idx,
                vcf_truth = truth_vcf_final,
                vcf_truth_idx = truth_vcf_final_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}",
                min_svlen_eval = min_svlen_eval,
                min_svlen_truth = min_svlen_truth,
                utils_docker = utils_docker,
                runtime_attr_subset_eval = runtime_attr_truvari_subset_eval,
                runtime_attr_subset_truth = runtime_attr_truvari_subset_truth,
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
                vcf_eval_idx = TruvariMatch.unmatched_vcf_idx,
                vcf_sv_truth = sv_truth_vcf_final,
                vcf_sv_truth_idx = sv_truth_vcf_final_idx,
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
                vcf_eval = eval_vcf_final,
                vcf_eval_idx = eval_vcf_final_idx,
                vcf_truth_snv = truth_vcf_final,
                vcf_truth_snv_idx = truth_vcf_final_idx,
                vcf_truth_sv = sv_truth_vcf_final,
                vcf_truth_sv_idx = sv_truth_vcf_final_idx,
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
                eval_vcf = eval_vcf_final,
                eval_vcf_idx = eval_vcf_final_idx,
                annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                matched_with_info_tsv = CollectMatchedIDsAndINFO.matched_with_info_tsv,
                eval_vep_header = ExtractEvalVepHeader.vep_header_txt,
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
        File benchmark_annotations_plots_tarball = MergePlotTarballs.merged_tarball
        File benchmark_annotations_summary_tsv = MergeBenchmarkSummaries.concatenated_tsv
        File benchmark_annotations_stats_tsv = MergeSummaryStats.concatenated_tsv
    }
}

task ExactMatch {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools isec \
            -c none \
            -n=2 \
            -p isec_matched \
            ~{vcf_eval} \
            ~{vcf_truth}
        
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            isec_matched/0000.vcf > eval_matched.tsv
        
        bcftools query \
            -f '%ID\n' \
            isec_matched/0001.vcf > truth_matched.tsv
        
        paste eval_matched.tsv truth_matched.tsv \
            | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"EXACT",$6,"SNV_indel"}' \
            > ~{prefix}.exact_matched.tsv
        
        bcftools isec \
            -C \
            -c none \
            -p isec_unmatched \
            ~{vcf_eval} \
            ~{vcf_truth}
        
        bgzip -c isec_unmatched/0000.vcf > ~{prefix}.unmatched.vcf.gz
        
        tabix -p vcf ~{prefix}.unmatched.vcf.gz
    >>>

    output {
        File annotation_tsv = "~{prefix}.exact_matched.tsv"
        File unmatched_vcf = "~{prefix}.unmatched.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}.unmatched.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 5,
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
        File vcf_eval_idx
        File vcf_truth_snv
        File vcf_truth_snv_idx
        File vcf_truth_sv
        File vcf_truth_sv_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess

annotation_tsv = "~{annotation_tsv}"
vcf_eval = "~{vcf_eval}"
vcf_truth_snv = "~{vcf_truth_snv}"
vcf_truth_sv = "~{vcf_truth_sv}"
prefix = "~{prefix}"

eval_to_truth = {}
eval_ids = set()
truth_ids = set()

with open(annotation_tsv) as f:
    for line in f:
        fields = line.strip().split('\t')
        eval_id = fields[4]
        truth_id = fields[6]
        eval_to_truth[eval_id] = truth_id
        eval_ids.add(eval_id)
        truth_ids.add(truth_id)

eval_info = {}
cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf_eval}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    parts = line.strip().split('\t', 1)
    if len(parts) == 2 and parts[0] in eval_ids:
        eval_info[parts[0]] = parts[1]
proc.wait()

truth_info = {}
for vcf in [vcf_truth_snv, vcf_truth_sv]:
    cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2 and parts[0] in truth_ids:
            truth_info[parts[0]] = parts[1]
    proc.wait()

with open(f"{prefix}.matched_with_info.tsv", 'w') as out:
    for eval_id, truth_id in eval_to_truth.items():
        eval_inf = eval_info.get(eval_id, '.')
        truth_inf = truth_info.get(truth_id, '.')
        out.write(f"{eval_id}\t{truth_id}\t{eval_inf}\t{truth_inf}\n")

EOF
    >>>

    output {
        File matched_with_info_tsv = "~{prefix}.matched_with_info.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) + 10,
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
        File vcf_idx
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
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
        
        if [ ! -s ~{matched_with_info_tsv} ] || [ $(cat ~{matched_with_info_tsv} | wc -l) -eq 0 ]; then
            echo "No matched variants found, creating empty shard"
            touch shards/matched.000000.tsv
        else
            cat ~{matched_with_info_tsv} | awk 'BEGIN{c=0;f=0} {print > sprintf("shards/matched.%06d.tsv", int(c/~{variants_per_shard})) ; c++} END{ }'
        fi
    >>>

    output {
        Array[File] shard_tsvs = glob("shards/*.tsv")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_with_info_tsv, "GB")) + 5,
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

        if [ ! -s ~{matched_shard_tsv} ] || [ $(cat ~{matched_shard_tsv} | wc -l) -eq 0 ]; then
            echo "Empty shard, creating empty output files"
            echo -e "af_key\teval_af\ttruth_af" > ~{prefix}.shard_~{shard_label}.af_pairs.tsv
            echo -e "category\teval\ttruth\tcount" > ~{prefix}.shard_~{shard_label}.vep_pairs.tsv
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
        File af_pairs_tsv = "~{prefix}.shard_~{shard_label}.af_pairs.tsv"
        File vep_pairs_tsv = "~{prefix}.shard_~{shard_label}.vep_pairs.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_shard_tsv, "GB")) + 5,
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
        mem_gb: 4,
        disk_gb: 2 * ceil(size(af_pair_tsvs[0], "GB")) + 5,
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
        File eval_vcf
        File eval_vcf_idx
        File annotation_tsv
        File matched_with_info_tsv
        File eval_vep_header
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
            --eval_vcf ~{eval_vcf} \
            --annotation_tsv ~{annotation_tsv} \
            --matched_with_info_tsv ~{matched_with_info_tsv} \
            --eval_vep_header ~{eval_vep_header} \
            --truth_vep_header ~{truth_vep_header}
    >>>

    output {
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(eval_vcf, "GB")) + 5,
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
        disk_gb: 2 * ceil(size(tarballs, "GB")) + 5,
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
