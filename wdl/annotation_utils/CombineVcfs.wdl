version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CombineVcfs {
    input {
        File a_vcf
        File a_vcf_idx
        File b_vcf
        File b_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_a
        RuntimeAttr? runtime_attr_subset_b
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_concat_merged
        RuntimeAttr? runtime_attr_concat_concatenated
    }

    call Helpers.CheckSampleMatch {
        input:
            vcf_a = a_vcf,
            vcf_a_idx = a_vcf_idx,
            vcf_b = b_vcf,
            vcf_b_idx = b_vcf_idx,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetA {
            input:
                vcf = a_vcf,
                vcf_idx = a_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.a.~{contig}.a",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_a
        }

        call Helpers.SubsetVcfToContig as SubsetB {
            input:
                vcf = b_vcf,
                vcf_idx = b_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.b",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_b
        }

        call Helpers.MergeVcfs {
            input:
                vcfs = [SubsetA.subset_vcf, SubsetB.subset_vcf],
                vcf_idxs = [SubsetA.subset_vcf_idx, SubsetB.subset_vcf_idx],
                prefix = "~{prefix}.~{contig}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call Helpers.ConcatVcfs {
            input:
                vcfs = [SubsetA.subset_vcf, SubsetB.subset_vcf],
                vcf_idxs = [SubsetA.subset_vcf_idx, SubsetB.subset_vcf_idx],
                prefix = "~{prefix}.~{contig}.concatenated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    call Helpers.ConcatVcfs as ConcatMerged {
        input:
            vcfs = MergeVcfs.merged_vcf,
            vcf_idxs = MergeVcfs.merged_vcf_idx,
            prefix = "~{prefix}.merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_merged
    }

    call Helpers.ConcatVcfs as ConcatConcatenated {
        input:
            vcfs = ConcatVcfs.concat_vcf,
            vcf_idxs = ConcatVcfs.concat_vcf_idx,
            prefix = "~{prefix}.concatenated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_concatenated
    }

    output {
        File merged_vcf = ConcatMerged.concat_vcf
        File merged_vcf_idx = ConcatMerged.concat_vcf_idx
        File concat_vcf = ConcatConcatenated.concat_vcf
        File concat_vcf_idx = ConcatConcatenated.concat_vcf_idx
    }
}
