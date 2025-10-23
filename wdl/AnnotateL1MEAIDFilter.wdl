version 1.0

import "general/Structs.wdl"

workflow AnnotateL1MEAIDFilter {
    input {
        File fa

        String prefix
        File rm_out

        String l1meaid_docker
        String l1meaid_filter_docker

        RuntimeAttr? runtime_attr_limeaid
        RuntimeAttr? runtime_attr_filter
    }

    call L1MEAID {
        input:
            fa = fa,
            prefix = prefix,
            rm_out = rm_out,
            docker = l1meaid_docker,
            runtime_attr_override = runtime_attr_limeaid
    }

    call L1MEAIDFilter {
        input:
            limeaid_output = L1MEAID.limeaid_output,
            prefix = prefix,
            docker = l1meaid_filter_docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File limeaid_output = L1MEAID.limeaid_output
        File l1meaid_filtered_output = L1MEAIDFilter.filtered_output
    }
}

task L1MEAID {
    input {
        File fa
        String prefix
        File rm_out
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/src/L1ME-AID/limeaid.py \
            -i ~{fa} \
            -r ~{rm_out} \
            -o ~{prefix}_limeaid.tsv
    >>>

    output {
        File limeaid_output = '~{prefix}_limeaid.tsv'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2*ceil(size(fa, "GB")) + 2*ceil(size(rm_out, "GB")) + 10,
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

task L1MEAIDFilter {
    input {
        File limeaid_output
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        perl /opt/src/utility/limeaid.filter.pl \
            ~{limeaid_output} \
            > ~{prefix}_filtered.tsv
    >>>

    output {
        File filtered_output = '~{prefix}_filtered.tsv'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2*ceil(size(limeaid_output, "GB")) + 5,
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
