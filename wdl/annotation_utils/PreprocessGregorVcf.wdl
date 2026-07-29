version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PreprocessGregorVcf {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File ref_fa
        File ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_shard_vcf
        RuntimeAttr? runtime_attr_normalize_vcf
        RuntimeAttr? runtime_attr_annotate_attributes
        RuntimeAttr? runtime_attr_rename_ids
        RuntimeAttr? runtime_attr_concat_shards
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard_vcf
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call Helpers.NormalizeVcf {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.normalized",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_vcf
            }

            call Helpers.AnnotateVariantAttributes {
                input:
                    vcf = NormalizeVcf.normalized_vcf,
                    vcf_idx = NormalizeVcf.normalized_vcf_idx,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate_attributes
            }

            call Helpers.RenameVariantIds {
                input:
                    vcf = AnnotateVariantAttributes.annotated_vcf,
                    vcf_idx = AnnotateVariantAttributes.annotated_vcf_idx,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.renamed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_ids
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = RenameVariantIds.renamed_vcf,
                    vcf_idxs = RenameVariantIds.renamed_vcf_idx,
                    allow_overlaps = true,
                    naive = false,
                    prefix = "~{prefix}.~{contig}.concatenated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_contig_vcf = select_first([ConcatShards.concat_vcf, RenameVariantIds.renamed_vcf[0]])
        File final_contig_vcf_idx = select_first([ConcatShards.concat_vcf_idx, RenameVariantIds.renamed_vcf_idx[0]])
    }

    output {
        Array[File] processed_vcfs = final_contig_vcf
        Array[File] processed_vcf_idxs = final_contig_vcf_idx
    }
}
