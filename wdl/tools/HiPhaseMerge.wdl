version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow HiPhaseMerge {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        Array[String] contigs
        String prefix

        String merge_args = "--merge id"

        String utils_docker

        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.MergeVcfs {
            input:
                vcfs = phased_vcfs,
                vcf_idxs = phased_vcf_idxs,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                extra_args = merge_args,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeVcfs.merged_vcf,
            vcfs_idx = MergeVcfs.merged_vcf_idx,
            prefix = "~{prefix}.hiphase_merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File hiphase_merged_vcf = ConcatVcfs.concat_vcf
        File hiphase_merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}