version 1.0

import "general/Structs.wdl"

workflow PALMERToVcf {
    input {
        Array[File] PALMER_calls
        Array[String] ME_type
        String sample
        File ref_fai

        String pipeline_docker
        RuntimeAttr? runtime_attr_override_palmer_to_vcf
        RuntimeAttr? runtime_attr_override_concat_sort_vcfs
    }

    scatter(i in range(length(PALMER_calls))) {
        call ConvertPALMERToVcf {
            input:
                PALMER_calls = PALMER_calls[i],
                ME_type       = ME_type[i],
                ref_fai       = ref_fai,
                sample        = sample,
                docker        = pipeline_docker,
                runtime_attr_override = runtime_attr_override_palmer_to_vcf
        }
    }

    call ConcatSortVcfs {
        input:
            vcfs   = ConvertPALMERToVcf.vcf,
            vcf_idxs = ConvertPALMERToVcf.vcf_idx,
            docker = pipeline_docker,
            sample = sample,
            runtime_attr_override = runtime_attr_override_concat_sort_vcfs
    }

    output {
        File PALMER_combined_vcf = ConcatSortVcfs.vcf
        File PALMER_combined_vcf_idx = ConcatSortVcfs.vcf_idx
    }
}

task ConvertPALMERToVcf {
    input {
        File PALMER_calls
        String ME_type
        String sample
        File ref_fai

        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(PALMER_calls, "GB") + size(ref_fai, "GB")) + 20

    command <<<
        set -euo pipefail

        python /opt/gnomad-lr/scripts/palmer/PALMER_to_VCF.py \
            ~{PALMER_calls} \
            ~{ME_type} \
            ~{sample} \
            ~{ref_fai} \
            | bcftools sort -Oz \
            > ~{sample}.PALMER_calls.~{ME_type}.vcf.gz
        
        tabix ~{sample}.PALMER_calls.~{ME_type}.vcf.gz
    >>>

    output {
        File vcf = "~{sample}.PALMER_calls.~{ME_type}.vcf.gz"
        File vcf_idx = "~{sample}.PALMER_calls.~{ME_type}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConcatSortVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String docker
        String sample
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 5*ceil(size(vcfs, "GB")) + 5

    command <<<
        set -euo pipefail

        bcftools concat -a ~{sep=" " vcfs} | bcftools sort -Oz > ~{sample}.PALMER_calls.merged.vcf.gz
        tabix ~{sample}.PALMER_calls.merged.vcf.gz
    >>>

    output {
        File vcf    = "~{sample}.PALMER_calls.merged.vcf.gz"
        File vcf_idx = "~{sample}.PALMER_calls.merged.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
