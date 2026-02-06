version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "TRGTMerge.wdl"

workflow HiPhaseMerge {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        Array[String] contigs
        String prefix

        String merge_args = "--merge id"

        File ref_fa
        File ref_fai

        String utils_docker
        String trgt_docker

        RuntimeAttr? runtime_attr_subset_trgt
        RuntimeAttr? runtime_attr_subset_integrated
        RuntimeAttr? runtime_attr_merge_trgt
        RuntimeAttr? runtime_attr_merge_integrated
        RuntimeAttr? runtime_attr_concat_trgt
        RuntimeAttr? runtime_attr_concat_integrated
    }

    scatter (i in range(length(phased_vcfs))) {
        call Helpers.SubsetVcfByArgs as SubsetWithTRID {
            input:
                vcf = phased_vcfs[i],
                vcf_idx = phased_vcf_idxs[i],
                include_args = "INFO/TRID != '.'",
                prefix = "~{prefix}.~{i}.trgt",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_trgt
        }

        call Helpers.SubsetVcfByArgs as SubsetWithoutTRID {
            input:
                vcf = phased_vcfs[i],
                vcf_idx = phased_vcf_idxs[i],
                exclude_args = "INFO/TRID != '.'",
                prefix = "~{prefix}.~{i}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_integrated
        }
    }

    scatter (contig in contigs) {
        call TRGTMerge.TRGTMergeContig as MergeTRGTVcfs {
            input:
                vcfs = SubsetWithTRID.subset_vcf,
                vcf_idxs = SubsetWithTRID.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.trgt",
                contig = contig,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                docker = trgt_docker,
                runtime_attr_override = runtime_attr_merge_trgt
        }

        call Helpers.MergeVcfs as MergeIntegratedVcfs {
            input:
                vcfs = SubsetWithoutTRID.subset_vcf,
                vcf_idxs = SubsetWithoutTRID.subset_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.integrated",
                extra_args = merge_args,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_integrated
        }
    }

    call Helpers.ConcatVcfs as ConcatTRGTVcfs {
        input:
            vcfs = MergeTRGTVcfs.merged_vcf,
            vcf_idxs = MergeTRGTVcfs.merged_vcf_idx,
            prefix = "~{prefix}.trgt",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_trgt
    }

    call Helpers.ConcatVcfs as ConcatIntegratedVcfs {
        input:
            vcfs = MergeIntegratedVcfs.merged_vcf,
            vcf_idxs = MergeIntegratedVcfs.merged_vcf_idx,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_integrated
    }

    output {
        File hiphase_merged_trgt_vcf = ConcatTRGTVcfs.concat_vcf
        File hiphase_merged_trgt_vcf_idx = ConcatTRGTVcfs.concat_vcf_idx
        File hiphase_merged_integrated_vcf = ConcatIntegratedVcfs.concat_vcf
        File hiphase_merged_integrated_vcf_idx = ConcatIntegratedVcfs.concat_vcf_idx
    }
}
