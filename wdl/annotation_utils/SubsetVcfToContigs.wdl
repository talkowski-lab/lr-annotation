version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow SubsetVcfToContigs {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = SubsetVcfToContig.subset_vcf,
            vcfs_idx = SubsetVcfToContig.subset_vcf_idx,
            prefix = "~{prefix}.subset_contigs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File subset_contigs_vcf = ConcatVcfs.concat_vcf
        File subset_contigs_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
