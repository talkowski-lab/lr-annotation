version 1.0

import "../utils/Structs.wdl"

workflow MethylationProfiling {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        String prefix

        String cpg_docker

        RuntimeAttr? runtime_attr_cpg_pileup
    }

    call CpgPileup {
        input:
            bam = bam,
            bai = bai,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = cpg_docker,
            runtime_attr_override = runtime_attr_cpg_pileup
    }

    output {
        File cpg_combined_bed = CpgPileup.pileup_combined_bed
        File cpg_combined_bed_idx = CpgPileup.pileup_combined_bed_idx
        File cpg_hap1_bed = CpgPileup.pileup_hap1_bed
        File cpg_hap1_bed_idx = CpgPileup.pileup_hap1_bed_idx
        File cpg_hap2_bed = CpgPileup.pileup_hap2_bed
        File cpg_hap2_bed_idx = CpgPileup.pileup_hap2_bed_idx
    }
}

task CpgPileup {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        aligned_bam_to_cpg_scores \
            --bam ~{bam} \
            --pileup-mode model \
            --ref ~{ref_fa} \
            --output-prefix ~{prefix} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --min-coverage 4 \
            --min-mapq 10
    >>>

    output {
        File pileup_combined_bed = "~{prefix}.combined.bed.gz"
        File pileup_combined_bed_idx = "~{prefix}.combined.bed.gz.tbi"
        File pileup_hap1_bed = "~{prefix}.hap1.bed.gz"
        File pileup_hap1_bed_idx = "~{prefix}.hap1.bed.gz.tbi"
        File pileup_hap2_bed = "~{prefix}.hap2.bed.gz"
        File pileup_hap2_bed_idx = "~{prefix}.hap2.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 12,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
