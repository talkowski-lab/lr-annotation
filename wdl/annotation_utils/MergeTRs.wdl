version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_subset_contig_tr
        RuntimeAttr? runtime_attr_integrate_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetBase {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.no_tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        call Helpers.SubsetVcfToContig as SubsetTR {
            input:
                vcf = tr_vcf,
                vcf_idx = tr_vcf,
                contig = contig,
                prefix = "~{prefix}.~{contig}.tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_tr
        }

        call Helpers.ConcatVcfs as IntegrateVcf {
            input:
                vcfs = [SubsetBase.subset_vcf, SubsetTR.subset_vcf],
                vcf_idxs = [SubsetBase.subset_vcf_idx, SubsetTR.subset_vcf_idx],
                allow_overlaps = true,
                naive = false,
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    call Helpers.ConcatVcfs as ConcatContigs {
        input:
            vcfs = IntegrateVcfs.concat_vcf,
            vcf_idxs = IntegrateVcfs.concat_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File merged_vcf = ConcatContigs.concat_vcf
        File merged_vcf_idx = ConcatContigs.concat_vcf_idx
    }
}
