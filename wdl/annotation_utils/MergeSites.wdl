version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeSites {
    input {
        File vcf
        File vcf_idx
        String prefix

        Float del_size_similarity = 0.8
        Float del_reciprocal_overlap = 0.8
        Int del_breakpoint_distance = 500
        Float del_sequence_similarity = 0.5
        Float del_sample_similarity = 0.5
        Int del_size_max = -1
        Int del_size_min = 0

        Float ins_reciprocal_overlap = 0.0
        Float ins_size_similarity = 0.5
        Int ins_breakpoint_distance = 10
        Float ins_sequence_similarity = 0.5
        Float ins_sample_similarity = 0.5
        Int ins_size_max = -1
        Int ins_size_min = 0

        String utils_docker

        RuntimeAttr? runtime_attr_subset_by_type
        RuntimeAttr? runtime_attr_collapse_dels
        RuntimeAttr? runtime_attr_collapse_ins
        RuntimeAttr? runtime_attr_concat
    }

    # Split VCF into DEL, INS and OTHER
    call Helpers.SubsetVcfByArgs as SubsetDels {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type~"del"',
            prefix = "~{prefix}.dels",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    call Helpers.SubsetVcfByArgs as SubsetIns {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type~"ins" || INFO/allele_type~"dup"',
            prefix = "~{prefix}.ins",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    call Helpers.SubsetVcfByArgs as SubsetOther {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type!~"del" && INFO/allele_type!~"ins" && INFO/allele_type!~"dup"',
            prefix = "~{prefix}.other",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    # Collapse DELs
    call Helpers.ConsolidateCollapsedSites as CollapseDels {
        input:
            vcf = SubsetDels.subset_vcf,
            vcf_idx = SubsetDels.subset_vcf_idx,
            pctovl = del_reciprocal_overlap,
            pctseq = del_sequence_similarity,
            pctsize = del_size_similarity,
            refdist = del_breakpoint_distance,
            sizemin = del_size_min,
            sizemax = del_size_max,
            keep_strategy = "common",
            sample_similarity = del_sample_similarity,
            set_merge_annotations = false,
            strip_format_to_gt = false,
            prefix = "~{prefix}.dels.collapsed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_collapse_dels
    }

    # Collapse INS
    call Helpers.ConsolidateCollapsedSites as CollapseIns {
        input:
            vcf = SubsetIns.subset_vcf,
            vcf_idx = SubsetIns.subset_vcf_idx,
            pctovl = ins_reciprocal_overlap,
            pctseq = ins_sequence_similarity,
            pctsize = ins_size_similarity,
            refdist = ins_breakpoint_distance,
            sizemin = ins_size_min,
            sizemax = ins_size_max,
            keep_strategy = "common",
            sample_similarity = ins_sample_similarity,
            set_merge_annotations = false,
            strip_format_to_gt = false,
            prefix = "~{prefix}.ins.collapsed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_collapse_ins
    }

    # Concat collapsed DELs + INS + OTHER
    call Helpers.ConcatVcfs {
        input:
            vcfs = [CollapseDels.consolidated_vcf, CollapseIns.consolidated_vcf, SubsetOther.subset_vcf],
            vcf_idxs = [CollapseDels.consolidated_vcf_idx, CollapseIns.consolidated_vcf_idx, SubsetOther.subset_vcf_idx],
            allow_overlaps = true,
            naive = false,
            prefix = "~{prefix}.merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File merged_vcf = ConcatVcfs.concat_vcf
        File merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
