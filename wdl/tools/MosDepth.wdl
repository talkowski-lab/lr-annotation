version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow MosDepth {
    input {
        File bam
        File bai
        Array[String] contigs
        String prefix

        Boolean quantize_mode

        String utils_docker

        RuntimeAttr? runtime_attr_run_mosdepth
    }

    scatter (contig in contigs) {
        call RunMosDepth {
            input:
                bam = bam,
                bai = bai,
                contig = contig,
                prefix = "~{prefix}.~{contig}.coverage",
                quantize_mode = quantize_mode,
                runtime_attr_override = runtime_attr_run_mosdepth
        }
    }

    output {
        Array[File] mosdepth_dist = RunMosDepth.dist
        Array[File] mosdepth_summary = RunMosDepth.summary
        Array[File] mosdepth_per_base = RunMosDepth.per_base
        Array[File] mosdepth_per_base_csi = RunMosDepth.per_base_csi
        Array[File] mosdepth_quantized_bed = select_all(RunMosDepth.quantized_bed)
        Array[File] mosdepth_quantized_bed_csi = select_all(RunMosDepth.quantized_bed_csi)
    }
}

task RunMosDepth {
    input {
        File bam
        File bai
        String contig
        String prefix
        Boolean quantize_mode
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ "~{quantize_mode}" == "true" ]; then
            export MOSDEPTH_Q0=NO_COVERAGE
            export MOSDEPTH_Q1=LOW_COVERAGE
            export MOSDEPTH_Q2=CALLABLE
            export MOSDEPTH_Q3=HIGH_COVERAGE

            mosdepth \
                -t 4 \
                -c "~{contig}" \
                -Q 1 \
                -x \
                --quantize 0:1:5:150: \
                ~{prefix} \
                ~{bam}
        else
            mosdepth \
                -t 2 \
                -c "~{contig}" \
                -Q 1 \
                -x \
                ~{prefix} \
                ~{bam}
        fi
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File per_base = "~{prefix}.per-base.bed.gz"
        File per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File? quantized_bed = "~{prefix}.quantized.bed.gz"
        File? quantized_bed_csi = "~{prefix}.quantized.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
