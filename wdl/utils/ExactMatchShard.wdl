version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExactMatchShard {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String region
        String prefix

        Int min_sv_length_truvari_vcf
        Int min_sv_length_truvari_truth_vcf
        String length_field_vcf
        String source_tag_truth_snv_indel_vcf

        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf

        String utils_docker

        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
    }

    call Helpers.SubsetVcfToRegion as SubsetExactEval {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            region = region,
            prefix = "~{prefix}.exact_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_exact_vcf
    }

    call Helpers.SubsetVcfToRegion as SubsetExactTruth {
        input:
            vcf = truth_snv_indel_vcf,
            vcf_idx = truth_snv_indel_vcf_idx,
            region = region,
            prefix = "~{prefix}.exact_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_exact_truth
    }

    if (defined(rename_id_string_vcf)) {
        call Helpers.RenameVariantIds as RenameEvalIds {
            input:
                vcf = SubsetExactEval.subset_vcf,
                vcf_idx = SubsetExactEval.subset_vcf_idx,
                prefix = "~{prefix}.exact_eval.renamed",
                id_format = select_first([rename_id_string_vcf]),
                strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_vcf
        }
    }

    if (defined(rename_id_string_truth_snv_indel_vcf)) {
        call Helpers.RenameVariantIds as RenameTruthIds {
            input:
                vcf = SubsetExactTruth.subset_vcf,
                vcf_idx = SubsetExactTruth.subset_vcf_idx,
                prefix = "~{prefix}.exact_truth.renamed",
                id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_truth
        }
    }

    File eval_vcf_final = select_first([RenameEvalIds.renamed_vcf, SubsetExactEval.subset_vcf])
    File eval_vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, SubsetExactEval.subset_vcf_idx])
    File truth_vcf_final = select_first([RenameTruthIds.renamed_vcf, SubsetExactTruth.subset_vcf])
    File truth_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, SubsetExactTruth.subset_vcf_idx])

    call Helpers.ExactMatch {
        input:
            vcf = eval_vcf_final,
            vcf_idx = eval_vcf_final_idx,
            truth_snv_indel_vcf = truth_vcf_final,
            truth_snv_indel_vcf_idx = truth_vcf_final_idx,
            source_tag = source_tag_truth_snv_indel_vcf,
            prefix = "~{prefix}.exact",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_exact_match
    }

    call Helpers.AppendAnnotationsFromVcf as AppendExactAnnotations {
        input:
            annotation_tsv = ExactMatch.annotation_tsv,
            truth_vcf = truth_vcf_final,
            truth_vcf_idx = truth_vcf_final_idx,
            is_sv_truth = false,
            prefix = "~{prefix}.exact_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_append_exact_annotations
    }

    call Helpers.SubsetVcfByLength as SubsetTruvariEval {
        input:
            vcf = ExactMatch.unmatched_vcf,
            vcf_idx = ExactMatch.unmatched_vcf_idx,
            length_field = length_field_vcf,
            min_length = min_sv_length_truvari_vcf,
            prefix = "~{prefix}.truvari_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_truvari_subset_vcf
    }

    call Helpers.SubsetVcfByArgs as SubsetTruvariTruth {
        input:
            vcf = truth_vcf_final,
            vcf_idx = truth_vcf_final_idx,
            include_args = "abs(ILEN) >= ~{min_sv_length_truvari_truth_vcf}",
            prefix = "~{prefix}.truvari_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_truvari_subset_truth
    }

    output {
        File annotated_tsv = AppendExactAnnotations.annotated_tsv
        File truvari_eval_vcf = SubsetTruvariEval.subset_vcf
        File truvari_eval_vcf_idx = SubsetTruvariEval.subset_vcf_idx
        File truvari_truth_vcf = SubsetTruvariTruth.subset_vcf
        File truvari_truth_vcf_idx = SubsetTruvariTruth.subset_vcf_idx
    }
}
