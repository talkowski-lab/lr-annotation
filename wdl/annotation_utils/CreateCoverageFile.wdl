version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow CreateCoverageFile {
    input {
        Array[Array[File]] mosdepth_bed_files
        Array[Array[File]] mosdepth_bed_indices
        Array[String] contigs
        String prefix

        Int bin_size
        Array[Int] thresholds

        RuntimeAttr? runtime_attr_bin
        RuntimeAttr? runtime_attr_merge
    }

    scatter (contig_idx in range(length(contigs))) {
        scatter (sample_idx in range(length(mosdepth_bed_files))) {
            File sample_contig_bed = mosdepth_bed_files[sample_idx][contig_idx]
            File sample_contig_bed_idx = mosdepth_bed_indices[sample_idx][contig_idx]
        }

        call ComputeBinnedCoverage {
            input:
                mosdepth_beds = sample_contig_bed,
                mosdepth_bed_indices = sample_contig_bed_idx,
                contig = contigs[contig_idx],
                bin_size = bin_size,
                thresholds = thresholds,
                prefix = "~{prefix}.~{contigs[contig_idx]}",
                runtime_attr_override = runtime_attr_bin
        }

        call FormatContigCoverage {
            input:
                binned_coverage = ComputeBinnedCoverage.binned_coverage,
                contig = contigs[contig_idx],
                thresholds = thresholds,
                prefix = "~{prefix}.~{contigs[contig_idx]}",
                runtime_attr_override = runtime_attr_bin
        }
    }

    call ConcatenateCoverages {
        input:
            formatted_coverages = FormatContigCoverage.formatted_coverage,
            thresholds = thresholds,
            prefix = "~{prefix}.coverage",
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File coverage = ConcatenateCoverages.merged_coverage
    }
}

task ComputeBinnedCoverage {
    input {
        Array[File] mosdepth_beds
        Array[File] mosdepth_bed_indices
        String contig
        Int bin_size
        Array[Int] thresholds
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
from statistics import mean, median
from collections import defaultdict

mosdepth_files = "~{sep=',' mosdepth_beds}".split(',')
bin_size = ~{bin_size}
thresholds = [~{sep=',' thresholds}]
contig = "~{contig}"
output_file = "~{prefix}.tsv"
num_samples = len(mosdepth_files)

bins = defaultdict(lambda: [0.0] * num_samples)

# Process each sample file
for sample_idx, bed_file in enumerate(mosdepth_files):
    # Track per-sample coverage for each bin
    sample_bin_coverages = defaultdict(list)
    
    with gzip.open(bed_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            coverage = float(parts[3])
            
            start_bin = (start // bin_size) * bin_size
            end_bin = ((end - 1) // bin_size) * bin_size
            
            for bin_start in range(start_bin, end_bin + 1, bin_size):
                bin_end = bin_start + bin_size
                overlap_start = max(start, bin_start)
                overlap_end = min(end, bin_end)
                overlap_length = overlap_end - overlap_start
                
                sample_bin_coverages[bin_start].extend([coverage] * overlap_length)
    
    # Compute per-sample mean for each bin and store
    for bin_start, coverages in sample_bin_coverages.items():
        bins[bin_start][sample_idx] = mean(coverages) if coverages else 0.0
    
    # Clear sample data to free memory
    sample_bin_coverages.clear()

# Write output
with open(output_file, 'w') as out:
    for bin_start in sorted(bins.keys()):
        bin_end = bin_start + bin_size
        sample_means = bins[bin_start]
        
        mean_cov = mean(sample_means)
        median_cov = int(median(sample_means))
        total_dp = sum(sample_means)
        over_counts = [sum(1 for sm in sample_means if sm > t) for t in thresholds]
        
        row = [contig, str(bin_start + 1), str(bin_end), str(mean_cov), str(median_cov), str(int(total_dp))]
        row.extend([str(c) for c in over_counts])
        out.write("\t".join(row) + "\n")
CODE

        bgzip ~{prefix}.tsv
    >>>

    output {
        File binned_coverage = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(mosdepth_beds, "GB")) + 20,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FormatContigCoverage {
    input {
        File binned_coverage
        String contig
        Array[Int] thresholds
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{binned_coverage} | awk -F'\t' 'BEGIN{OFS="\t"} {printf "%s:%s-%s", $1, $2, $3; for(i=4; i<=NF; i++) printf "\t%s", $i; printf "\n"}' | bgzip > ~{prefix}.formatted.tsv.gz
    >>>

    output {
        File formatted_coverage = "~{prefix}.formatted.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(binned_coverage, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatenateCoverages {
    input {
        Array[File] formatted_coverages
        Array[Int] thresholds
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{sep=' ' formatted_coverages} | bgzip > ~{prefix}.tsv.gz
    >>>

    output {
        File merged_coverage = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(formatted_coverages, "GB")) + 10,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
