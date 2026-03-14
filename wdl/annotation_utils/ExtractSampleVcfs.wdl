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

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_subset_snv
        RuntimeAttr? runtime_attr_subset_sv
        RuntimeAttr? runtime_attr_concat_snv
        RuntimeAttr? runtime_attr_concat_sv
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetContig {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        scatter (sample_id in sample_ids) {
            call Helpers.ExtractSample {
                input:
                    vcf = SubsetContig.subset_vcf,
                    vcf_idx = SubsetContig.subset_vcf_idx,
                    sample = sample_id,
                    prefix = "~{prefix}.~{sample_id}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_sample
            }

            call Helpers.SubsetVcfByArgs as SubsetSnvIndel {
                input:
                    vcf = ExtractSample.subset_vcf,
                    vcf_idx = ExtractSample.subset_vcf_idx,
                    include_args = 'abs(INFO/allele_length) < ~{min_sv_length} && INFO/SOURCE = \"DeepVariant\"',
                    prefix = "~{prefix}.~{sample_id}.~{contig}.snv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_snv
            }

            call Helpers.SubsetVcfByArgs as SubsetSv {
                input:
                    vcf = ExtractSample.subset_vcf,
                    vcf_idx = ExtractSample.subset_vcf_idx,
                    include_args = 'abs(INFO/allele_length) >= ~{min_sv_length} && INFO/SOURCE = \"HPRC_SV_Integration\"',
                    prefix = "~{prefix}.~{sample_id}.~{contig}.sv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_sv
            }
        }
    }

    Array[Array[File]] snv_indel_vcfs_by_sample = transpose(SubsetSnvIndel.subset_vcf)
    Array[Array[File]] snv_indel_vcf_idxs_by_sample = transpose(SubsetSnvIndel.subset_vcf_idx)
    Array[Array[File]] sv_vcfs_by_sample = transpose(SubsetSv.subset_vcf)
    Array[Array[File]] sv_vcf_idxs_by_sample = transpose(SubsetSv.subset_vcf_idx)

    scatter (i in range(length(sample_ids))) {
        call Helpers.ConcatVcfs as ConcatSampleSnv {
            input:
                vcfs = snv_indel_vcfs_by_sample[i],
                vcf_idxs = snv_indel_vcf_idxs_by_sample[i],
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{sample_ids[i]}.snv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_snv
        }

        call Helpers.ConcatVcfs as ConcatSampleSv {
            input:
                vcfs = sv_vcfs_by_sample[i],
                vcf_idxs = sv_vcf_idxs_by_sample[i],
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{sample_ids[i]}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_sv
        }
    }

    output {
        Array[File] snv_indel_vcfs = ConcatSampleSnv.concat_vcf
        Array[File] snv_indel_vcf_idxs = ConcatSampleSnv.concat_vcf_idx
        Array[File] sv_vcfs = ConcatSampleSv.concat_vcf
        Array[File] sv_vcf_idxs = ConcatSampleSv.concat_vcf_idx
    }
}
