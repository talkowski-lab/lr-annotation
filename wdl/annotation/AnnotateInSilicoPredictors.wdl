version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateInSilicoPredictors {

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht

        String genome_build = "GRCh38"
        String annotate_in_silico_predictors_script

        String hail_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call AnnotateInSilicoPredictorsTask {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.in_silico_predictors",
                cadd_ht = cadd_ht,
                pangolin_ht = pangolin_ht,
                phylop_ht = phylop_ht,
                revel_ht = revel_ht,
                spliceai_ht = spliceai_ht,
                docker = hail_docker,
                annotate_in_silico_predictors_script = annotate_in_silico_predictors_script,
                genome_build = genome_build,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateInSilicoPredictorsTask.annotated_vcf,
            vcf_idxs = AnnotateInSilicoPredictorsTask.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.in_silico_predictors",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File insilico_annotated_vcf = ConcatVcfs.concat_vcf
        File insilico_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateInSilicoPredictorsTask {
    input {
        File vcf
        File vcf_idx

        String prefix
        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht

        String docker
        String annotate_in_silico_predictors_script
        String genome_build = "GRCh38"

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(10.0 + size(vcf, "GB") * 10.0),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        curl ~{annotate_in_silico_predictors_script} > annotate.py

        python3 annotate.py \
            --build ~{genome_build} \
            --cadd_ht ~{cadd_ht} \
            --pangolin_ht ~{pangolin_ht} \
            --phylop_ht ~{phylop_ht} \
            --revel_ht ~{revel_ht} \
            --spliceai_ht ~{spliceai_ht} \
            --vcf ~{vcf} \
            --output_vcf ~{prefix}.vcf.gz
        
        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
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