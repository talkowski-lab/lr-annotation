version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int min_svlen
        File coding_gtf
        File noncoding_bed

        String utils_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_annotate_func
        RuntimeAttr? runtime_attr_concat_unannotated
        RuntimeAttr? runtime_attr_concat_annotated
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_postprocess
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfBySize as SubsetVcfAnnotated {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                locus = contig,
                min_size = min_svlen,
                prefix = "~{prefix}.~{contig}.annotate",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ConvertToSymbolic as PreprocessVcf {
            input:
                vcf = SubsetVcfAnnotated.subset_vcf,
                vcf_idx = SubsetVcfAnnotated.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.preprocessed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_preprocess
        }

        call AnnotateFunctionalConsequences {
            input:
                vcf = PreprocessVcf.processed_vcf,
                vcf_idx = PreprocessVcf.processed_vcf_idx,
                noncoding_bed = noncoding_bed,
                coding_gtf = coding_gtf,
                prefix = "~{prefix}.~{contig}.functionally_annotated",
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_annotate_func
        }

        call Helpers.RevertSymbolicAlleles as PostprocessVcf {
            input:
                annotated_vcf = AnnotateFunctionalConsequences.anno_vcf,
                annotated_vcf_idx = AnnotateFunctionalConsequences.anno_vcf_idx,
                original_vcf = SubsetVcfAnnotated.subset_vcf,
                original_vcf_idx = SubsetVcfAnnotated.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.postprocessed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_postprocess
        }

        call Helpers.ExtractVcfAnnotations as ExtractAnnotations {
            input:
                vcf = PostprocessVcf.reverted_vcf,
                vcf_idx = PostprocessVcf.reverted_vcf_idx,
                original_vcf = SubsetVcfAnnotated.subset_vcf,
                original_vcf_idx = SubsetVcfAnnotated.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker
        }
    }

    call Helpers.ConcatTsvs as MergeTsvs {
        input:
            tsvs = ExtractAnnotations.annotations_tsv,
            prefix = prefix + ".svannotate_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_annotated
    }

    call Helpers.MergeHeaderLines as MergeHeaders {
        input:
            header_files = ExtractAnnotations.annotations_header,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File annotations_tsv_svannotate = MergeTsvs.concatenated_tsv
        File annotations_header_svannotate = MergeHeaders.merged_header
    }
}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File noncoding_bed
        File coding_gtf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = 1000 * ceil(0.8 * select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
            -V ~{vcf} \
            --non-coding-bed ~{noncoding_bed} \
            --protein-coding-gtf ~{coding_gtf} \
            -O ~{prefix}.vcf.gz
    >>>

    output {
        File anno_vcf = "~{prefix}.vcf.gz"
        File anno_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
