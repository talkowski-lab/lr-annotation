version 1.0

import "Structs.wdl"

workflow AnnotateL1MEAID {
    input {
        File fasta
        File rm_out
        String l1meaid_docker
        String l1meaid_filter_docker
        
        RuntimeAttr? runtime_attr_override_limeaid
        RuntimeAttr? runtime_attr_override_filter
    }

    call L1MEAID {
        input:
            fasta = fasta,
            rm_out = rm_out,
            docker = l1me_aid_docker,
            runtime_attr_override = runtime_attr_override_limeaid
    }

    call L1MEAIDFilter {
        input:
            limeaid_output = L1MEAID.limeaid_output,
            docker = l1meaid_filter_docker,
            runtime_attr_override = runtime_attr_override_filter
    }

    output {
        File limeaid_output = L1MEAID.limeaid_output
        File filtered_output = L1MEAIDFilter.filtered_output
    }
}

task L1MEAID {
    input {
        File fasta
        File rm_out
        String docker
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fasta, "GB")) + 2*ceil(size(rm_out, "GB")) + 10
    String prefix = basename(fasta, ".fasta")

    command <<<
        set -euo pipefail

        python /opt/src/L1ME-AID/limeaid.py \
            -i ~{fasta} \
            -r ~{rm_out} \
            -o ~{prefix}_limeaid.tsv
    >>>

    output {
        File limeaid_output = '~{prefix}_limeaid.tsv'
    }

    RuntimeAttr runtime_default = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: select_first([runtime_attr.docker, runtime_default.docker])
    }
}

task L1MEAIDFilter {
    input {
        File limeaid_output
        String docker
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(limeaid_output, "GB")) + 5
    String prefix = basename(limeaid_output, ".tsv")

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
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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
