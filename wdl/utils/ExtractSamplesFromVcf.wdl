version 1.0

import "../utils/Structs.wdl"

workflow ExtractSamplesFromVcf {
    input {
        File vcf
        File vcf_idx
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_extract_samples
    }

    call ExtractSamples {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_extract_samples
    }

    output {
        File sample_list = ExtractSamples.sample_list
    }
}

task ExtractSamples {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools query -l ~{vcf} | grep -v '^$' > ~{prefix}.samples.txt
    >>>

    output {
        File sample_list = "~{prefix}.samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
