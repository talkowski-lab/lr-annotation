version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillPhasedGenotypes {
    input {
        File phased_vcf
        File phased_vcf_idx
        File unphased_vcf
        File unphased_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker
        
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_split
        RuntimeAttr? runtime_attr_fill
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetPhased {
            input:
                vcf = phased_vcf,
                vcf_idx = phased_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.phased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.SubsetVcfToContig as SubsetUnphased {
            input:
                vcf = unphased_vcf,
                vcf_idx = unphased_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.unphased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.FillGenotypes {
            input:
                phased_vcf = SubsetPhased.subset_vcf,
                phased_vcf_idx = SubsetPhased.subset_vcf_idx,
                unphased_vcf = SubsetUnphased.subset_vcf,
                unphased_vcf_idx = SubsetUnphased.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.filled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_fill
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = FillGenotypes.filled_vcf,
            vcf_idxs = FillGenotypes.filled_vcf_idx,
            prefix = "~{prefix}.filled",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File hiphase_phased_integrated_vcf = ConcatVcfs.concat_vcf
        File hiphase_phased_integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
