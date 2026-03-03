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
        File cpg_combined_bed = CpgPileup.pileup_bed_combined
        File cpg_combined_bw = CpgPileup.pileup_bigwig_combined
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

    RuntimeAttr default_attr = object {
        cpu_cores: 12,
        mem_gb: 48,
        disk_gb: ceil((size(bam, "GB") + size(ref_fa, "GB")) * 2 + 20),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    command <<<
        set -euo pipefail
        
        aligned_bam_to_cpg_scores \
            --threads ~{threads} \
            --bam ~{bam} \
            --ref ~{ref_fa} \
            --output-prefix ~{prefix} \
            --min-mapq 1 \
            --min-coverage 10 \
            --model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite

        gzip ~{prefix}.*.bed
    >>>

    output {
        File pileup_bed_combined = "~{prefix}.combined.bed.gz"
        File pileup_bigwig_combined = "~{prefix}.combined.bw"
    }
}
