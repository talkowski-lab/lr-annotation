version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow CreateCoverageFile {
    input {
        Array[File] mosdepth_bed_files
        Array[File] mosdepth_bed_indices
        Array[String] contigs
        String prefix

        Int bin_size
        Array[Int] thresholds

        String utils_docker

        RuntimeAttr? runtime_attr_compute_coverage
        RuntimeAttr? runtime_attr_merge_coverages
    }

    scatter (contig in contigs) {
        call ComputeBinnedCoverage {
            input:
                mosdepth_beds = mosdepth_bed_files,
                mosdepth_bed_indices = mosdepth_bed_indices,
                contig = contig,
                bin_size = bin_size,
                thresholds = thresholds,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_compute_coverage
        }
    }

    call ConcatenateCoverages {
        input:
            formatted_coverages = ComputeBinnedCoverage.binned_coverage,
            thresholds = thresholds,
            prefix = "~{prefix}.coverage",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_coverages
    }

    output {
        File binned_coverage_tsv = ConcatenateCoverages.merged_coverage
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
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
import subprocess
from statistics import mean, median

mosdepth_files = "~{sep=',' mosdepth_beds}".split(',')
bin_size = ~{bin_size}
thresholds = [~{sep=',' thresholds}]
contig = "~{contig}"
output_file = "~{prefix}.tsv"
num_samples = len(mosdepth_files)

iterators = []
current_rows = []

for bed_file in mosdepth_files:
    proc = subprocess.Popen(['tabix', '-0', bed_file, contig], stdout=subprocess.PIPE, text=True)
    iterators.append(proc.stdout)
    
    line = proc.stdout.readline()
    if line:
        parts = line.strip().split('\t')
        current_rows.append((int(parts[1]), int(parts[2]), float(parts[3])))
    else:
        current_rows.append(None)

min_pos = min((row[0] for row in current_rows if row), default=None)
with gzip.open(output_file + '.gz', 'wt') as out:
    current_bin = (min_pos // bin_size) * bin_size
    
    while any(row is not None for row in current_rows):
        bin_end = current_bin + bin_size
        bin_coverages = [[] for _ in range(num_samples)]
        has_data = False
        
        for sample_idx in range(num_samples):
            while current_rows[sample_idx] is not None:
                start, end, coverage = current_rows[sample_idx]
                
                if end <= current_bin:
                    # Row before bin, skip it
                    line = iterators[sample_idx].readline()
                    if line:
                        parts = line.strip().split('\t')
                        current_rows[sample_idx] = (int(parts[1]), int(parts[2]), float(parts[3]))
                    else:
                        current_rows[sample_idx] = None
                elif start >= bin_end:
                    # Row after bin, keep it for next bin
                    break
                else:
                    # Row overlaps bin
                    has_data = True
                    overlap_start = max(start, current_bin)
                    overlap_end = min(end, bin_end)
                    overlap_length = overlap_end - overlap_start
                    bin_coverages[sample_idx].extend([coverage] * overlap_length)
                    
                    if end > bin_end:
                        # Row extends beyond bin, keep for next
                        break
                    else:
                        # Row consumed, get next
                        line = iterators[sample_idx].readline()
                        if line:
                            parts = line.strip().split('\t')
                            current_rows[sample_idx] = (int(parts[1]), int(parts[2]), float(parts[3]))
                        else:
                            current_rows[sample_idx] = None
        
        if has_data:
            sample_means = [mean(cov) if cov else 0.0 for cov in bin_coverages]
            mean_cov = mean(sample_means)
            median_cov = int(median(sample_means))
            total_dp = sum(sample_means)
            over_counts = [sum(1 for sm in sample_means if sm > t) for t in thresholds]
            
            row = [contig, str(current_bin + 1), str(bin_end), str(mean_cov), str(median_cov), str(int(total_dp))]
            row.extend(map(str, over_counts))
            out.write("\t".join(row) + "\n")
        
        current_bin += bin_size
    
    for it in iterators:
        it.close()
CODE
    >>>

    output {
        File binned_coverage = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(mosdepth_beds, "GB")) + 20,
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

task ConcatenateCoverages {
    input {
        Array[File] formatted_coverages
        Array[Int] thresholds
        String prefix
        String docker
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
        mem_gb: 2,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
