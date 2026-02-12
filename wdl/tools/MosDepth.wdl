version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AlignedMetrics {
    parameter_meta {
        aligned_bam: "Aligned BAM file"
        aligned_bai: "Index for aligned BAM file"
        ref_fasta: "Reference FASTA file"
        ref_dict: "Reference dictionary file"
        gcs_output_dir: "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."
    }

    input {
        File aligned_bam
        File aligned_bai

        File ref_fasta
        File ref_dict

        String prefix

        String? gcs_output_dir

        RuntimeAttr? runtime_attr_read_metrics
        RuntimeAttr? runtime_attr_make_chr_interval_list
        RuntimeAttr? runtime_attr_mosdepth
        RuntimeAttr? runtime_attr_summarize_depth
        RuntimeAttr? runtime_attr_flag_stats
        RuntimeAttr? runtime_attr_finalize
    }

    call ReadMetrics as AlignedReadMetrics {
        input:
            bam = aligned_bam,
            prefix = prefix + ".read_metrics",
            runtime_attr_override = runtime_attr_read_metrics
    }

    call MakeChrIntervalList {
        input:
            ref_dict = ref_dict,
            prefix = prefix + ".chr_intervals",
            runtime_attr_override = runtime_attr_make_chr_interval_list
    }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call MosDepth {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                chr = chr_info[0],
                prefix = "~{prefix}.~{chr_info[0]}.coverage",
                runtime_attr_override = runtime_attr_mosdepth
        }

        call SummarizeDepth {
            input:
                regions = MosDepth.regions,
                prefix = "~{prefix}.~{chr_info[0]}.coverage_summary",
                runtime_attr_override = runtime_attr_summarize_depth
        }
    }

    call FlagStats as AlignedFlagStats {
        input:
            bam = aligned_bam,
            prefix = prefix + ".flag_stats",
            runtime_attr_override = runtime_attr_flag_stats
    }

    if (defined(gcs_output_dir)) {
        String outdir = sub(gcs_output_dir + "", "/$", "")

        call Helpers.FinalizeToDir as FFYieldAligned {
            input:
                outdir = outdir + "/yield_aligned/",
                files = [
                    AlignedFlagStats.flag_stats,
                    AlignedReadMetrics.np_hist,
                    AlignedReadMetrics.range_gap_hist,
                    AlignedReadMetrics.zmw_hist,
                    AlignedReadMetrics.prl_counts,
                    AlignedReadMetrics.prl_hist,
                    AlignedReadMetrics.prl_nx,
                    AlignedReadMetrics.prl_yield_hist,
                    AlignedReadMetrics.rl_counts,
                    AlignedReadMetrics.rl_hist,
                    AlignedReadMetrics.rl_nx,
                    AlignedReadMetrics.rl_yield_hist
                ],
                runtime_attr_override = runtime_attr_finalize
        }

        call Helpers.FinalizeToDir as FFCoverageFullDist {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.full_dist,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageGlobalDist {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.global_dist,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageRegionDist {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.region_dist,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageRegions {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.regions,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageRegionsCsi {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.regions_csi,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageQuantizedDist {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.quantized_dist,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageQuantized {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.quantized,
                runtime_attr_override = runtime_attr_finalize
        }
        call Helpers.FinalizeToDir as FFCoverageQuantizedCsi {
            input:
                outdir = outdir + "/coverage/",
                files = MosDepth.quantized_csi,
                runtime_attr_override = runtime_attr_finalize
        }

        call Helpers.FinalizeToDir as FFDepthSummaries {
            input:
                outdir = outdir + "/coverage_summaries/",
                files = SummarizeDepth.cov_summary,
                runtime_attr_override = runtime_attr_finalize
        }
    }

    output {
        File aligned_flag_stats = AlignedFlagStats.flag_stats
        Array[File] coverage_full_dist = MosDepth.full_dist
        Array[File] coverage_global_dist = MosDepth.global_dist
        Array[File] coverage_region_dist = MosDepth.region_dist
        Array[File] coverage_regions = MosDepth.regions
        Array[File] coverage_regions_csi = MosDepth.regions_csi
        Array[File] coverage_quantized_dist = MosDepth.quantized_dist
        Array[File] coverage_quantized = MosDepth.quantized
        Array[File] coverage_quantized_csi  = MosDepth.quantized_csi
        File aligned_np_hist = AlignedReadMetrics.np_hist
        File aligned_range_gap_hist = AlignedReadMetrics.range_gap_hist
        File aligned_zmw_hist = AlignedReadMetrics.zmw_hist
        File aligned_prl_counts = AlignedReadMetrics.prl_counts
        File aligned_prl_hist = AlignedReadMetrics.prl_hist
        File aligned_prl_nx = AlignedReadMetrics.prl_nx
        File aligned_prl_yield_hist = AlignedReadMetrics.prl_yield_hist
        File aligned_rl_counts = AlignedReadMetrics.rl_counts
        File aligned_rl_hist = AlignedReadMetrics.rl_hist
        File aligned_rl_nx = AlignedReadMetrics.rl_nx
        File aligned_rl_yield_hist = AlignedReadMetrics.rl_yield_hist
        File raw_chr_intervals = MakeChrIntervalList.raw_chrs
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
        mem_gb: 1,
        disk_gb: 10,
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

task MosDepth {
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

        mosdepth -t 4 \
            -c "~{chr}" \
            -Q 1 \
            -n \
            -x \
            ~{prefix}.full \
            ~{bam}
        
        mosdepth \
            -t 4 \
            -c "~{chr}" \
            -Q 1 \
            -b ~{window_size} \
            -n \
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
            -n \
            -x \
            --quantize 0:1:5:150: \
            ~{prefix}.quantized \
            ~{bam}
    >>>

    output {
        File full_dist = "~{prefix}.full.mosdepth.global.dist.txt"
        File global_dist = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist = "~{prefix}.mosdepth.region.dist.txt"
        File regions = "~{prefix}.regions.bed.gz"
        File regions_csi = "~{prefix}.regions.bed.gz.csi"
        File quantized_dist = "~{prefix}.quantized.mosdepth.global.dist.txt"
        File quantized = "~{prefix}.quantized.quantized.bed.gz"
        File quantized_csi = "~{prefix}.quantized.quantized.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(bam, "GB") + size(bai, "GB")) + 10,
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
        mem_gb: 2,
        disk_gb: 2 * ceil(size(regions, "GB")),
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
        disk_gb: ceil(size(bam, "GB")),
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

task CallableLoci {
    parameter_meta {
        bam_file: "Input BAM file to analyze"
        bam_index: "Index file for input BAM" 
        ref_fasta: "Reference FASTA file"
        ref_fasta_index: "Index file for reference FASTA"
        ref_dict: "Dictionary file for reference FASTA"
        prefix: "Prefix for output files"
        min_depth: "Minimum depth for a locus to be considered callable"
        min_base_quality: "Minimum base quality for a base to be considered"
        min_mapping_quality: "Minimum mapping quality for a read to be considered"
        runtime_attr_override: "Runtime attributes override struct"
    }

    input {
        File bam_file
        File bam_index
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String prefix
        Int min_depth = 5
        Int min_base_quality = 20
        Int min_mapping_quality = 10
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gatk CallableLoci \
            -R ~{ref_fasta} \
            -I ~{bam_file} \
            -O ~{prefix}.callable_status.bed \
            --summary ~{prefix}.callable_status.summary.txt \
            --min-depth ~{min_depth} \
            --min-base-quality ~{min_base_quality} \
            --min-mapping-quality ~{min_mapping_quality}
    >>>

    output {
        File callable_loci_bed = "~{prefix}.callable_status.bed"
        File callable_loci_summary = "~{prefix}.callable_status.summary.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(bam_file, "GiB") + size(ref_fasta, "GiB") + 20),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:jonn-fa5c895f58b729f00589f5deb23d56efb929aa1d-4.6.1.0-5-gfa5c895f5"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
