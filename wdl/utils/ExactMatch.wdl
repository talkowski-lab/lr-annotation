version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExactMatch {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String contig
        String prefix

        Int? shard_bin_size_exact_match

        Int min_sv_length_truvari_vcf
        Int min_sv_length_truvari_truth_vcf
        String length_field_vcf
        String source_tag_truth_snv_indel_vcf

        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf

        String utils_docker

        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_concat_truvari_eval
        RuntimeAttr? runtime_attr_concat_truvari_truth
    }

    if (defined(shard_bin_size_exact_match)) {
        call Helpers.CreateContigShards as CreateExactShards {
            input:
                vcfs = [vcf, truth_snv_indel_vcf],
                vcf_idxs = [vcf_idx, truth_snv_indel_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size_exact_match]),
                prefix = "~{prefix}.exact_shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_exact_shards
        }

        scatter (k in range(length(CreateExactShards.shard_regions))) {
            call Helpers.SubsetVcfToRegion as SubsetExactEvalRaw {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    region = CreateExactShards.shard_regions[k],
                    prefix = "~{prefix}.exact_eval_~{k}.raw",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_exact_vcf
            }

            call Helpers.StripGenotypes as StripGenotypesShard {
                input:
                    vcf = SubsetExactEvalRaw.subset_vcf,
                    vcf_idx = SubsetExactEvalRaw.subset_vcf_idx,
                    prefix = "~{prefix}.exact_eval_~{k}.stripped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_strip_genotypes
            }

            call Helpers.SubsetVcfToRegion as SubsetExactTruth {
                input:
                    vcf = truth_snv_indel_vcf,
                    vcf_idx = truth_snv_indel_vcf_idx,
                    region = CreateExactShards.shard_regions[k],
                    prefix = "~{prefix}.exact_truth_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_exact_truth
            }

            if (defined(rename_id_string_vcf)) {
                call Helpers.RenameVariantIds as RenameEvalIdsShard {
                    input:
                        vcf = StripGenotypesShard.stripped_vcf,
                        vcf_idx = StripGenotypesShard.stripped_vcf_idx,
                        prefix = "~{prefix}.exact_eval_~{k}.renamed",
                        id_format = select_first([rename_id_string_vcf]),
                        strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_rename_vcf
                }
            }

            if (defined(rename_id_string_truth_snv_indel_vcf)) {
                call Helpers.RenameVariantIds as RenameTruthIdsShard {
                    input:
                        vcf = SubsetExactTruth.subset_vcf,
                        vcf_idx = SubsetExactTruth.subset_vcf_idx,
                        prefix = "~{prefix}.exact_truth_~{k}.renamed",
                        id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                        strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_rename_truth
                }
            }

            File exact_eval_vcf_final = select_first([RenameEvalIdsShard.renamed_vcf, StripGenotypesShard.stripped_vcf])
            File exact_eval_vcf_final_idx = select_first([RenameEvalIdsShard.renamed_vcf_idx, StripGenotypesShard.stripped_vcf_idx])
            File exact_truth_vcf_final = select_first([RenameTruthIdsShard.renamed_vcf, SubsetExactTruth.subset_vcf])
            File exact_truth_vcf_final_idx = select_first([RenameTruthIdsShard.renamed_vcf_idx, SubsetExactTruth.subset_vcf_idx])

            call Helpers.ExactMatch as ExactMatchShard {
                input:
                    vcf = exact_eval_vcf_final,
                    vcf_idx = exact_eval_vcf_final_idx,
                    truth_snv_indel_vcf = exact_truth_vcf_final,
                    truth_snv_indel_vcf_idx = exact_truth_vcf_final_idx,
                    source_tag = source_tag_truth_snv_indel_vcf,
                    prefix = "~{prefix}.exact_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_exact_match
            }

            call Helpers.AppendAnnotationsFromVcf as AppendExactAnnotationsShard {
                input:
                    annotation_tsv = ExactMatchShard.annotation_tsv,
                    truth_vcf = exact_truth_vcf_final,
                    truth_vcf_idx = exact_truth_vcf_final_idx,
                    is_sv_truth = false,
                    prefix = "~{prefix}.exact_annotated_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_append_exact_annotations
            }

            call Helpers.SubsetVcfByLength as SubsetTruvariEvalShard {
                input:
                    vcf = ExactMatchShard.unmatched_vcf,
                    vcf_idx = ExactMatchShard.unmatched_vcf_idx,
                    length_field = length_field_vcf,
                    min_length = min_sv_length_truvari_vcf,
                    prefix = "~{prefix}.truvari_eval_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_truvari_subset_vcf
            }

            call Helpers.SubsetVcfByArgs as SubsetTruvariTruthShard {
                input:
                    vcf = exact_truth_vcf_final,
                    vcf_idx = exact_truth_vcf_final_idx,
                    include_args = "abs(ILEN) >= ~{min_sv_length_truvari_truth_vcf}",
                    prefix = "~{prefix}.truvari_truth_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_truvari_subset_truth
            }
        }

        call Helpers.ConcatTsvs as ConcatExactAnnotations {
            input:
                tsvs = AppendExactAnnotationsShard.annotated_tsv,
                sort_output = true,
                preserve_header = true,
                prefix = "~{prefix}.exact_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_exact_annotations
        }

        call Helpers.ConcatVcfs as ConcatTruvariEval {
            input:
                vcfs = SubsetTruvariEvalShard.subset_vcf,
                vcf_idxs = SubsetTruvariEvalShard.subset_vcf_idx,
                allow_overlaps = false,
                naive = false,
                prefix = "~{prefix}.truvari_eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_truvari_eval
        }

        call Helpers.ConcatVcfs as ConcatTruvariTruth {
            input:
                vcfs = SubsetTruvariTruthShard.subset_vcf,
                vcf_idxs = SubsetTruvariTruthShard.subset_vcf_idx,
                allow_overlaps = false,
                naive = false,
                prefix = "~{prefix}.truvari_truth",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_truvari_truth
        }
    }

    if (!defined(shard_bin_size_exact_match)) {
        call Helpers.StripGenotypes as StripGenotypesFull {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                prefix = "~{prefix}.eval.stripped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_strip_genotypes
        }

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIdsFull {
                input:
                    vcf = StripGenotypesFull.stripped_vcf,
                    vcf_idx = StripGenotypesFull.stripped_vcf_idx,
                    prefix = "~{prefix}.eval.renamed",
                    id_format = select_first([rename_id_string_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_vcf
            }
        }

        if (defined(rename_id_string_truth_snv_indel_vcf)) {
            call Helpers.RenameVariantIds as RenameTruthIdsFull {
                input:
                    vcf = truth_snv_indel_vcf,
                    vcf_idx = truth_snv_indel_vcf_idx,
                    prefix = "~{prefix}.truth.renamed",
                    id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_truth
            }
        }

        File vcf_final = select_first([RenameEvalIdsFull.renamed_vcf, StripGenotypesFull.stripped_vcf])
        File vcf_final_idx = select_first([RenameEvalIdsFull.renamed_vcf_idx, StripGenotypesFull.stripped_vcf_idx])
        File truth_snv_indel_vcf_final = select_first([RenameTruthIdsFull.renamed_vcf, truth_snv_indel_vcf])
        File truth_snv_indel_vcf_final_idx = select_first([RenameTruthIdsFull.renamed_vcf_idx, truth_snv_indel_vcf_idx])

        call Helpers.ExactMatch as ExactMatchFull {
            input:
                vcf = vcf_final,
                vcf_idx = vcf_final_idx,
                truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                source_tag = source_tag_truth_snv_indel_vcf,
                prefix = "~{prefix}.exact",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call Helpers.AppendAnnotationsFromVcf as AppendExactAnnotationsFull {
            input:
                annotation_tsv = ExactMatchFull.annotation_tsv,
                truth_vcf = ExactMatchFull.matched_truth_vcf,
                truth_vcf_idx = ExactMatchFull.matched_truth_vcf_idx,
                is_sv_truth = false,
                prefix = "~{prefix}.exact_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_append_exact_annotations
        }

        call Helpers.SubsetVcfByLength as SubsetTruvariEvalFull {
            input:
                vcf = ExactMatchFull.unmatched_vcf,
                vcf_idx = ExactMatchFull.unmatched_vcf_idx,
                length_field = length_field_vcf,
                min_length = min_sv_length_truvari_vcf,
                prefix = "~{prefix}.truvari_eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_subset_vcf
        }

        call Helpers.SubsetVcfByArgs as SubsetTruvariTruthFull {
            input:
                vcf = truth_snv_indel_vcf_final,
                vcf_idx = truth_snv_indel_vcf_final_idx,
                include_args = "abs(ILEN) >= ~{min_sv_length_truvari_truth_vcf}",
                prefix = "~{prefix}.truvari_truth",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_subset_truth
        }
    }

    output {
        File annotated_tsv = select_first([ConcatExactAnnotations.concatenated_tsv, AppendExactAnnotationsFull.annotated_tsv])
        File truvari_eval_vcf = select_first([ConcatTruvariEval.concat_vcf, SubsetTruvariEvalFull.subset_vcf])
        File truvari_eval_vcf_idx = select_first([ConcatTruvariEval.concat_vcf_idx, SubsetTruvariEvalFull.subset_vcf_idx])
        File truvari_truth_vcf = select_first([ConcatTruvariTruth.concat_vcf, SubsetTruvariTruthFull.subset_vcf])
        File truvari_truth_vcf_idx = select_first([ConcatTruvariTruth.concat_vcf_idx, SubsetTruvariTruthFull.subset_vcf_idx])
    }
}
