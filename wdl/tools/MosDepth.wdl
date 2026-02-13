version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow MosDepth {
    input {
        File bam
        File bai
        String prefix

        File ref_fasta
        File ref_dict

        RuntimeAttr? runtime_attr_read_metrics
        RuntimeAttr? runtime_attr_make_chr_interval_list
        RuntimeAttr? runtime_attr_run_mosdepth
        RuntimeAttr? runtime_attr_summarize_depth
        RuntimeAttr? runtime_attr_flag_stats
    }

    call ReadMetrics {
        input:
            bam = bam,
            prefix = "~{prefix}.read_metrics",
            runtime_attr_override = runtime_attr_read_metrics
    }

    call MakeChrIntervalList {
        input:
            ref_dict = ref_dict,
            prefix = "~{prefix}.intervals",
            runtime_attr_override = runtime_attr_make_chr_interval_list
    }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call RunMosDepth {
            input:
                bam = bam,
                bai = bai,
                chr = chr_info[0],
                prefix = "~{prefix}.~{chr_info[0]}.coverage",
                runtime_attr_override = runtime_attr_run_mosdepth
        }

        call SummarizeDepth {
            input:
                regions = RunMosDepth.regions,
                prefix = "~{prefix}.~{chr_info[0]}.summary",
                runtime_attr_override = runtime_attr_summarize_depth
        }
    }

    call FlagStats {
        input:
            bam = bam,
            prefix = "~{prefix}.flag_stats",
            runtime_attr_override = runtime_attr_flag_stats
    }

    output {
        File mosdepth_np_hist = ReadMetrics.np_hist
        File mosdepth_range_gap_hist = ReadMetrics.range_gap_hist
        File mosdepth_zmw_hist = ReadMetrics.zmw_hist
        File mosdepth_prl_counts = ReadMetrics.prl_counts
        File mosdepth_prl_hist = ReadMetrics.prl_hist
        File mosdepth_prl_nx = ReadMetrics.prl_nx
        File mosdepth_prl_yield_hist = ReadMetrics.prl_yield_hist
        File mosdepth_rl_counts = ReadMetrics.rl_counts
        File mosdepth_rl_hist = ReadMetrics.rl_hist
        File mosdepth_rl_nx = ReadMetrics.rl_nx
        File mosdepth_rl_yield_hist = ReadMetrics.rl_yield_hist
        File mosdepth_raw_chr_intervals = MakeChrIntervalList.raw_chrs

        Array[File] mosdepth_dist = RunMosDepth.dist
        Array[File] mosdepth_summary = RunMosDepth.summary
        Array[File] mosdepth_per_base = RunMosDepth.per_base
        Array[File] mosdepth_per_base_csi = RunMosDepth.per_base_csi
        Array[File] mosdepth_regions = RunMosDepth.regions
        Array[File] mosdepth_regions_csi = RunMosDepth.regions_csi
        Array[File] mosdepth_quantized_dist = RunMosDepth.quantized_dist
        Array[File] mosdepth_quantized_summary = RunMosDepth.quantized_summary
        Array[File] mosdepth_quantized_per_base  = RunMosDepth.quantized_per_base
        Array[File] mosdepth_quantized_per_base_csi = RunMosDepth.quantized_per_base_csi
        Array[File] mosdepth_quantized_regions = RunMosDepth.quantized_regions
        Array[File] mosdepth_quantized_regions_csi = RunMosDepth.quantized_regions_csi
        Array[File] mosdepth_quantized_bed = RunMosDepth.quantized_bed
        Array[File] mosdepth_quantized_bed_csi = RunMosDepth.quantized_bed_csi

        Array[File] mosdepth_summarized_coverage = SummarizeDepth.cov_summary
        File mosdepth_flag_stats = FlagStats.flag_stats
    }
}

task ReadMetrics {
    input {
        File bam
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        java -jar /usr/local/bin/gatk.jar ComputeLongReadMetrics \
            -I ~{bam} \
            -O ~{prefix} \
            -DF WellformedReadFilter
    >>>

    output {
        File np_hist = "~{prefix}.np_hist.txt"
        File range_gap_hist = "~{prefix}.range_gap_hist.txt"
        File zmw_hist = "~{prefix}.zmw_hist.txt"
        File prl_counts = "~{prefix}.prl_counts.txt"
        File prl_hist = "~{prefix}.prl_hist.txt"
        File prl_nx = "~{prefix}.prl_nx.txt"
        File prl_yield_hist = "~{prefix}.prl_yield_hist.txt"
        File rl_counts = "~{prefix}.rl_counts.txt"
        File rl_hist = "~{prefix}.rl_hist.txt"
        File rl_nx = "~{prefix}.rl_nx.txt"
        File rl_yield_hist = "~{prefix}.rl_yield_hist.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 50,
        disk_gb: 2 * ceil(size(bam, "GB")),
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MakeChrIntervalList {
    input {
        File ref_dict
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        grep '^@SQ' ~{ref_dict} \
            | awk '{ print $2 "\t" 1 "\t" $3 }' \
            | sed 's/[SL]N://g' \
            | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV \
            > ~{prefix}.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("~{prefix}.txt")
        File raw_chrs = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(ref_dict, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task RunMosDepth {
    input {
        File bam
        File bai
        String chr
        String prefix
        Int window_size = 500
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mosdepth \
            -t 4 \
            -c "~{chr}" \
            -Q 1 \
            -x \
            ~{prefix} \
            ~{bam}

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth \
            -t 4 \
            -c "~{chr}" \
            -Q 1 \
            -x \
            --quantize 0:1:5:150: \
            ~{prefix}.quantized \
            ~{bam}
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File per_base = "~{prefix}.per-base.bed.gz"
        File per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File regions = "~{prefix}.regions.bed.gz"
        File regions_csi = "~{prefix}.regions.bed.gz.csi"
        File quantized_dist = "~{prefix}.quantized.mosdepth.global.dist.txt"
        File quantized_summary = "~{prefix}.quantized.mosdepth.summary.txt"
        File quantized_per_base = "~{prefix}.quantized.per-base.bed.gz"
        File quantized_per_base_csi = "~{prefix}.quantized.per-base.bed.gz.csi"
        File quantized_regions = "~{prefix}.quantized.regions.bed.gz"
        File quantized_regions_csi = "~{prefix}.quantized.regions.bed.gz.csi"
        File quantized_bed = "~{prefix}.quantized.bed.gz"
        File quantized_bed_csi = "~{prefix}.quantized.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(bam, "GB")) + 10,
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

task SummarizeDepth {
    input {
        File regions
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ((echo 'chr start stop cov_mean cov_sd cov_q1 cov_median cov_q3 cov_iqr') && \
         (zcat ~{regions} | datamash first 1 first 2 last 3 mean 4 sstdev 4 q1 4 median 4 q3 4 iqr 4)) | \
         column -t > ~{prefix}.txt
    >>>

    output {
        File cov_summary = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(regions, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FlagStats {
    input {
        File bam
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools flagstat \
            ~{bam} \
            > ~{prefix}.txt
    >>>

    output {
        File flag_stats = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bam, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
