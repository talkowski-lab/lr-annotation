version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow MosDepth {
    input {
        File bam
        File bai
        Array[String] contigs
        String prefix

        Int? bin_size

        File? ref_fa
        File? ref_fai

        String mosdepth_docker

        RuntimeAttr? runtime_attr_run_mosdepth
    }

    scatter (contig in contigs) {
        call RunMosDepth {
            input:
                bam = bam,
                bai = bai,
                contig = contig,
                bin_size = bin_size,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.coverage",
                docker = mosdepth_docker,
                runtime_attr_override = runtime_attr_run_mosdepth
        }
    }

    output {
        Array[File] mosdepth_dist = RunMosDepth.dist
        Array[File] mosdepth_summary = RunMosDepth.summary
        Array[File] mosdepth_per_base = select_all(RunMosDepth.per_base)
        Array[File] mosdepth_per_base_csi = select_all(RunMosDepth.per_base_csi)
        Array[File] mosdepth_regions_bed = select_all(RunMosDepth.regions_bed)
        Array[File] mosdepth_regions_bed_csi = select_all(RunMosDepth.regions_bed_csi)
    }
}

task RunMosDepth {
    input {
        File bam
        File bai
        String contig
        Int? bin_size
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mosdepth \
            -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -c "~{contig}" \
            -x \
            ~{if defined(ref_fa) then "-f " + ref_fa else ""} \
            ~{if defined(bin_size) then "--by " + bin_size + " --no-per-base" else ""} \
            ~{prefix} \
            ~{bam}
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File? per_base = "~{prefix}.per-base.bed.gz"
        File? per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File? regions_bed = "~{prefix}.regions.bed.gz"
        File? regions_bed_csi = "~{prefix}.regions.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + 10,
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
