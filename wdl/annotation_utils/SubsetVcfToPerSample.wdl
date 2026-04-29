version 1.0

import "../utils/Helpers.wdl"

workflow SubsetVcfToPerSample {
    input {
        Array[File] cohort_vcfs
        Array[File] cohort_vcf_idxs
        Array[String] sample_ids
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (i in range(length(contigs))) {
        String contig = contigs[i]

        scatter (sample_id in sample_ids) {
            call Helpers.ExtractSample {
                input:
                    vcf = cohort_vcfs[i],
                    vcf_idx = cohort_vcf_idxs[i],
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
                allow_overlaps = true,
                naive = false,
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
