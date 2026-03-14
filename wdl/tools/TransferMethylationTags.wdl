version 1.0

import "../utils/Structs.wdl"

workflow TransferMethylationTags {
    input {
        File unaligned_bam
        File aligned_bam
        File aligned_bai
        String prefix
        
        String docker 

        RuntimeAttr? runtime_attr_transfer
    }

    call TransferTags {
        input:
            unaligned_bam = unaligned_bam,
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            prefix = "~{prefix}.tagged",
            docker = docker,
            runtime_attr_override = runtime_attr_transfer
    }

    output {
        File methylation_tagged_bam = TransferTags.tagged_bam
        File methylation_tagged_bai = TransferTags.tagged_bai
    }
}

task TransferTags {
    input {
        File unaligned_bam
        File aligned_bam
        File aligned_bai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "Extracting read names from unaligned BAM..."
        samtools view ~{unaligned_bam} \
            | cut -f1 \
            | sort -u \
            > subset_qnames.txt

        echo "Subsetting aligned BAM using extracted read names..."
        samtools view \
            -h \
            -N \
            subset_qnames.txt \
            ~{aligned_bam} \
            -o subset_aligned.bam

        echo "Transferring MM and ML tags..."
        python3 <<CODE
import array
import pysam

tags_dict = {}

# Read the unaligned subset BAM and cache the tags in memory
with pysam.AlignmentFile("~{unaligned_bam}", "rb", check_sq=False) as ubam:
    for read in ubam:
        try:
            mm = read.get_tag('MM')
            ml = list(read.get_tag('ML'))
            tags_dict[read.query_name] = (mm, ml)
        except KeyError:
            pass # Skip reads that don't have methylation tags

# Iterate over the subset aligned BAM, add tags and write out
with pysam.AlignmentFile("subset_aligned.bam", "rb") as abam:
    with pysam.AlignmentFile("~{prefix}.unsorted.bam", "wb", header=abam.header) as outbam:
        for read in abam:
            # if read.is_supplementary or read.is_secondary:
            #     continue
            if read.query_name in tags_dict:
                mm, ml = tags_dict[read.query_name]
                read.set_tag('MM', mm)
                read.set_tag('ML', array.array('B', ml))
                outbam.write(read)
CODE

        echo "Sorting and indexing..."
        samtools sort -@ 4 -o "~{prefix}.bam" "~{prefix}.unsorted.bam"
        samtools index -@ 4 "~{prefix}.bam"
    >>>

    output {
        File tagged_bam = "~{prefix}.bam"
        File tagged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(unaligned_bam, "GB") + size(aligned_bam, "GB")) + 20,
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
