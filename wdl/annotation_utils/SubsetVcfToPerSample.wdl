version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetVcfToPerSample {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
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
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    sample = sample_id,
                    normalize_output = true,
                    prefix = "~{prefix}.~{sample_id}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_sample
            }
        }
    }

    # Transpose from [contig][sample] to [sample][contig] for per-sample concat
    Array[Array[File]] vcfs_by_sample = transpose(ExtractSample.subset_vcf)
    Array[Array[File]] vcf_idxs_by_sample = transpose(ExtractSample.subset_vcf_idx)

    scatter (i in range(length(sample_ids))) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = vcfs_by_sample[i],
                vcf_idxs = vcf_idxs_by_sample[i],
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{sample_ids[i]}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        Array[File] subset_vcfs = ConcatVcfs.concat_vcf
        Array[File] subset_vcf_idxs = ConcatVcfs.concat_vcf_idx
    }
}
