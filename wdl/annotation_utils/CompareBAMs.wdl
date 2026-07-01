version 1.0

import "../utils/Structs.wdl"

workflow CompareBAMs {
    input {
        File bam1
        File bam2
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_compare_bams
    }

    call CompareBamReadCounts {
        input:
            bam1 = bam1,
            bam2 = bam2,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare_bams
    }

    output {
        File comparison_tsv = CompareBamReadCounts.comparison_tsv
    }
}

task CompareBamReadCounts {
    input {
        File bam1
        File bam2
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{bam1} \
            | awk 'BEGIN{OFS="\t"} {print $1, ($10 != "*") ? length($10) : 0}' \
            | LC_ALL=C sort -k1,1 -T . -S 4G > bam1_reads.tsv &
        PID1=$!

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{bam2} \
            | awk 'BEGIN{OFS="\t"} {print $1, ($10 != "*") ? length($10) : 0}' \
            | LC_ALL=C sort -k1,1 -T . -S 4G > bam2_reads.tsv &
        PID2=$!

        wait $PID1
        wait $PID2

        BAM1_COUNT=$(wc -l < bam1_reads.tsv)
        BAM2_COUNT=$(wc -l < bam2_reads.tsv)

        # Streaming merge-join on read name: outputs "name\tlen1\tlen2" for matched IDs
        LC_ALL=C join -t $'\t' -1 1 -2 1 bam1_reads.tsv bam2_reads.tsv > matched_reads.tsv

        MATCHED_COUNT=$(wc -l < matched_reads.tsv)
        SAME_LEN_COUNT=$(awk -F'\t' '$2 == $3 {n++} END {print n+0}' matched_reads.tsv)

        printf "metric\tvalue\n" > ~{prefix}.bam_comparison.tsv
        printf "bam1_total_reads\t%s\n" "$BAM1_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "bam2_total_reads\t%s\n" "$BAM2_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "matched_reads\t%s\n" "$MATCHED_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "matched_reads_with_same_sequence_length\t%s\n" "$SAME_LEN_COUNT" >> ~{prefix}.bam_comparison.tsv
    >>>

    output {
        File comparison_tsv = "~{prefix}.bam_comparison.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: ceil(size(bam1, "GB") + size(bam2, "GB")) + 20,
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
