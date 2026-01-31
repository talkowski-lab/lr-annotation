version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExtractSampleVcfs {
    input {
        Array[String] sample_ids
        File cohort_vcf
        File cohort_vcf_idx
        String prefix

        Int min_sv_length
        String? extra_args

        String utils_docker

        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_subset_snv
        RuntimeAttr? runtime_attr_subset_sv
    }

    scatter (sample_id in sample_ids) {
        call Helpers.ExtractSample {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                sample = sample_id,
                extra_args = extra_args,
                prefix = "~{prefix}.~{sample_id}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_sample
        }

        call Helpers.SubsetVcfByLength as SubsetSnv {
            input:
                vcf = ExtractSample.subset_vcf,
                vcf_idx = ExtractSample.subset_vcf_idx,
                max_size = min_sv_length - 1,
                prefix = "~{prefix}.~{sample_id}.snv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_snv
        }

        call Helpers.SubsetVcfByLength as SubsetSv {
            input:
                vcf = ExtractSample.subset_vcf,
                vcf_idx = ExtractSample.subset_vcf_idx,
                min_length = min_sv_length,
                prefix = "~{prefix}.~{sample_id}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_sv
        }
    }

    output {
        Array[File] snv_vcfs = SubsetSnv.subset_vcf
        Array[File] snv_vcf_idxs = SubsetSnv.subset_vcf_idx
        Array[File] sv_vcfs = SubsetSv.subset_vcf
        Array[File] sv_vcf_idxs = SubsetSv.subset_vcf_idx
    }
}
