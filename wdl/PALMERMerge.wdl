version 1.0

import "general/Structs.wdl"

workflow PALMERMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs

        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_override_merge_vcfs
    }

    call MergeVCFs {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override_merge_vcfs
    }

    output {
        File palmer_merged_vcf = MergeVCFs.merged_vcf
        File palmer_merged_vcf_idx = MergeVCFs.merged_vcf_idx
    }
}

task MergeVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir inputs

        vcf_files=(~{sep=" " vcfs})
        vcf_idxs=(~{sep=" " vcf_idxs})

        for i in ${!vcf_files[@]}; do
            ln -s ${vcf_files[$i]} inputs/input_$i.vcf.gz
            ln -s ${vcf_idxs[$i]} inputs/input_$i.vcf.gz.tbi
        done

        bcftools merge \
            --missing-to-ref \
            -Oz -o tmp.merged.vcf.gz \
            inputs/*.vcf.gz

        bcftools annotate \
            --set-id '%INFO/ME_TYPE\_%CHROM\_%POS\_%INFO/SVLEN' \
            -Oz -o ~{prefix}.merged.vcf.gz \
            tmp.merged.vcf.gz

        tabix ~{prefix}.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_idx = "~{prefix}.merged.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2*ceil(size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
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
