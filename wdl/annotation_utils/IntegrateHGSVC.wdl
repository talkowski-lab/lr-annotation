version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateHGSVC {
    input {
        File snv_vcf
        File snv_vcf_idx
        File indel_vcf
        File indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        String snv_source_tag
        String snv_source_tag_description
        String indel_source_tag
        String indel_source_tag_description
        String sv_source_tag
        String sv_source_tag_description
        
        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_add_info
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = [snv_vcf, indel_vcf, sv_vcf],
            vcfs_idx = [snv_vcf_idx, indel_vcf_idx, sv_vcf_idx],
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        # SNV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSnv {
            input:
                vcf = snv_vcf,
                vcf_idx = snv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSnv {
            input:
                vcf = SubsetSnv.subset_vcf,
                vcf_idx = SubsetSnv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoSnv {
            input:
                vcf = AnnotateSnv.annotated_vcf,
                vcf_idx = AnnotateSnv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = snv_source_tag,
                tag_description = snv_source_tag_description,
                prefix = "~{prefix}.~{contig}.snv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # Indel VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetIndel {
            input:
                vcf = indel_vcf,
                vcf_idx = indel_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateIndel {
            input:
                vcf = SubsetIndel.subset_vcf,
                vcf_idx = SubsetIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoIndel {
            input:
                vcf = AnnotateIndel.annotated_vcf,
                vcf_idx = AnnotateIndel.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = indel_source_tag,
                tag_description = indel_source_tag_description,
                prefix = "~{prefix}.~{contig}.indel.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # SV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSv {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_idx = AnnotateSv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = sv_source_tag,
                tag_description = sv_source_tag_description,
                prefix = "~{prefix}.~{contig}.sv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [AddInfoSnv.annotated_vcf, AddInfoIndel.annotated_vcf, AddInfoSv.annotated_vcf],
                vcfs_idx = [AddInfoSnv.annotated_vcf_idx, AddInfoIndel.annotated_vcf_idx, AddInfoSv.annotated_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }
    
    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcfs_idx = MergeContigVcfs.concat_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}