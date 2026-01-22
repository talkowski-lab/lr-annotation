version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow SubsetVcfToSamples {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Array[String] samples

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToSampleList {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                samples = samples,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = SubsetVcfToSampleList.subset_vcf,
            vcfs_idx = SubsetVcfToSampleList.subset_vcf_idx,
            merge_sort = true,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File subset_samples_vcf = ConcatVcfs.concat_vcf
        File subset_samples_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
