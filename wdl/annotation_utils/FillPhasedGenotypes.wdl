version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

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
                prefix = prefix + ".phased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.SubsetVcfToContig as SubsetUnphased {
            input:
                vcf = unphased_vcf,
                vcf_idx = unphased_vcf_idx,
                contig = contig,
                prefix = prefix + ".unphased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.SplitMultiallelics {
            input:
                vcf = SubsetPhased.subset_vcf,
                vcf_idx = SubsetPhased.subset_vcf_idx,
                prefix = prefix + ".phased.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_split
        }

        call Helpers.FillGenotypesFromUnphased {
            input:
                phased_vcf = SplitMultiallelics.split_vcf,
                phased_vcf_idx = SplitMultiallelics.split_vcf_idx,
                unphased_vcf = SubsetUnphased.subset_vcf,
                unphased_vcf_idx = SubsetUnphased.subset_vcf_idx,
                prefix = prefix + "." + contig + ".filled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_fill
        }

        call Helpers.AnnotateVariantAttributes {
            input:
                vcf = FillGenotypesFromUnphased.filled_vcf,
                vcf_idx = FillGenotypesFromUnphased.filled_vcf_idx,
                prefix = prefix + "." + contig + ".filled.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateVariantAttributes.annotated_vcf,
            vcfs_idx = AnnotateVariantAttributes.annotated_vcf_idx,
            prefix = prefix + ".filled.annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File filled_vcf = ConcatVcfs.concat_vcf
        File filled_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}