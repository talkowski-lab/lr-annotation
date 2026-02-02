version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateVcfs {
    input {
        File snv_indel_vcf
        File snv_indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        Int min_sv_length = 50
        String snv_indel_vcf_source_tag
        String snv_indel_vcf_size_flag
        String snv_indel_vcf_size_flag_description
        String sv_vcf_source_tag
        String sv_vcf_size_flag
        String sv_vcf_size_flag_description
        File? swap_samples_snv_indel
        File? swap_samples_sv
        
        String utils_docker
        
        RuntimeAttr? runtime_attr_swap_samples_snv_indel
        RuntimeAttr? runtime_attr_swap_samples_sv
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_filter_snv_indel
        RuntimeAttr? runtime_attr_filter_sv
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_annotate_attributes
    }

    if (defined(swap_samples_snv_indel)) {
        call Helpers.SwapSampleIds as SwapSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                sample_swap_list = select_first([swap_samples_snv_indel]),
                prefix = prefix + ".snv_indel.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples_snv_indel
        }
    }

    if (defined(swap_samples_sv)) {
        call Helpers.SwapSampleIds as SwapSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                sample_swap_list = select_first([swap_samples_sv]),
                prefix = prefix + ".sv.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples_sv
        }
    }
    
    File final_snv_indel_vcf = select_first([SwapSnvIndel.swapped_vcf, snv_indel_vcf])
    File final_snv_indel_vcf_idx = select_first([SwapSnvIndel.swapped_vcf_idx, snv_indel_vcf_idx])
    File final_sv_vcf = select_first([SwapSv.swapped_vcf, sv_vcf])
    File final_sv_vcf_idx = select_first([SwapSv.swapped_vcf_idx, sv_vcf_idx])

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = [final_snv_indel_vcf, final_sv_vcf],
            vcfs_idx = [final_snv_indel_vcf, final_sv_vcf],
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        # SNV Indel VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSnvIndel {
            input:
                vcf = final_snv_indel_vcf,
                vcf_idx = final_snv_indel_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv_indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.SplitMultiallelics as SplitSnvIndel {
            input:
                vcf = SubsetSnvIndel.subset_vcf,
                vcf_idx = SubsetSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSnvIndel {
            input:
                vcf = SplitSnvIndel.split_vcf,
                vcf_idx = SplitSnvIndel.split_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_attributes
        }

        call Helpers.AddInfo as AddInfoSnvIndel {
            input:
                vcf = AnnotateSnvIndel.annotated_vcf,
                vcf_idx = AnnotateSnvIndel.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = snv_indel_vcf_source_tag,
                tag_description = "Source of variant call",
                prefix = "~{prefix}.~{contig}.snv_indel.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.AddFilter as AddFilterSnvIndel {
            input:
                vcf = AddInfoSnvIndel.annotated_vcf,
                vcf_idx = AddInfoSnvIndel.annotated_vcf_idx,
                filter_name = snv_indel_vcf_size_flag,
                filter_description = snv_indel_vcf_size_flag_description,
                filter_expression = "abs(INFO/allele_length) >= ~{min_sv_length}",
                prefix = "~{prefix}.~{contig}.snv_indel.flagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        # SV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = final_sv_vcf,
                vcf_idx = final_sv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.SplitMultiallelics as SplitSv {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSv {
            input:
                vcf = SplitSv.split_vcf,
                vcf_idx = SplitSv.split_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_attributes
        }

        call Helpers.AddInfo as AddInfoSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_idx = AnnotateSv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = sv_vcf_source_tag,
                tag_description = "Source of variant call",
                prefix = "~{prefix}.~{contig}.sv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.AddFilter as AddFilterSv {
            input:
                vcf = AddInfoSv.annotated_vcf,
                vcf_idx = AddInfoSv.annotated_vcf_idx,
                filter_name = sv_vcf_size_flag,
                filter_description = sv_vcf_size_flag_description,
                filter_expression = "abs(INFO/allele_length) < ~{min_sv_length}",
                prefix = "~{prefix}.~{contig}.sv.flagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [AddFilterSnvIndel.flagged_vcf, AddFilterSv.flagged_vcf],
                vcfs_idx = [AddFilterSnvIndel.flagged_vcf_idx, AddFilterSv.flagged_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcfs_idx = MergeContigVcfs.concat_vcf_idx,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
