version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExtractSampleVcfs {
    input {
        Array[String] sample_ids
        Array[String] contigs
        File cohort_vcf
        File cohort_vcf_idx
        String prefix

        Int min_sv_length

        String utils_docker

        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_subset_snv
        RuntimeAttr? runtime_attr_subset_sv
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (sample_id in sample_ids) {
        scatter (contig in contigs) {
            call Helpers.ExtractSample as ExtractSampleContig {
                input:
                    vcf = cohort_vcf,
                    vcf_idx = cohort_vcf_idx,
                    sample = sample_id,
                    extra_args = "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.~{sample_id}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_sample
            }

            call Helpers.SubsetVcfByArgs as SubsetSnvIndel {
                input:
                    vcf = ExtractSampleContig.subset_vcf,
                    vcf_idx = ExtractSampleContig.subset_vcf_idx,
                    include_args = 'abs(INFO/allele_length) < ~{min_sv_length} && INFO/SOURCE = \"DeepVariant\"',
                    prefix = "~{prefix}.~{contig}.~{sample_id}.snv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_snv
            }

            call Helpers.SubsetVcfByArgs as SubsetSv {
                input:
                    vcf = ExtractSampleContig.subset_vcf,
                    vcf_idx = ExtractSampleContig.subset_vcf_idx,
                    include_args = 'abs(INFO/allele_length) >= ~{min_sv_length} && INFO/SOURCE = \"HPRC_SV_Integration\"',
                    prefix = "~{prefix}.~{contig}.~{sample_id}.sv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_sv
            }
        }

        call Helpers.ConcatVcfs as ConcatSnvIndels {
            input:
                vcfs = SubsetSnvIndel.subset_vcf,
                vcf_idxs = SubsetSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{sample_id}.snv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }

        call Helpers.ConcatVcfs as ConcatSvs {
            input:
                vcfs = SubsetSv.subset_vcf,
                vcf_idxs = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{sample_id}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        Array[File] snv_indel_vcfs = ConcatSnvIndels.concat_vcf
        Array[File] snv_indel_vcf_idxs = ConcatSnvIndels.concat_vcf_idx
        Array[File] sv_vcfs = ConcatSvs.concat_vcf
        Array[File] sv_vcf_idxs = ConcatSvs.concat_vcf_idx
    }
}
