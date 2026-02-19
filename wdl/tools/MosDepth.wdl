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
        RuntimeAttr? runtime_attr_concat_tsvs
        RuntimeAttr? runtime_attr_compress_index
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

    call Helpers.ConcatTsvs as ConcatPerBase {
        input:
            tsvs = RunMosDepth.per_base,
            prefix = "~{prefix}.per-base",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    call CompressAndIndex {
        input:
            bed_file = ConcatPerBase.concatenated_tsv,
            prefix = "~{prefix}.per-base",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compress_index
    }

    output {
        Array[File] mosdepth_dist = RunMosDepth.dist
        Array[File] mosdepth_summary = RunMosDepth.summary
        File mosdepth_per_base = CompressAndIndex.compressed_bed
        File mosdepth_per_base_csi = CompressAndIndex.compressed_bed_csi
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
            export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
            export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
            export MOSDEPTH_Q2=CALLABLE      # 5..149
            export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

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
        mem_gb: 2,
        disk_gb: ceil(size(bam, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CompressAndIndex {
    input {
        File bed_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bgzip -c ~{bed_file} > ~{prefix}.bed.gz

        tabix -p bed -C ~{prefix}.bed.gz
    >>>

    output {
        File compressed_bed = "~{prefix}.bed.gz"
        File compressed_bed_csi = "~{prefix}.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_file, "GB")) + 10,
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
