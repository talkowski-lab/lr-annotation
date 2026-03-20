version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TransferMethylationTags {
    input {
        Array[String] unaligned_bam_paths
        File aligned_bam
        File aligned_bai
        Array[String] contigs
        String prefix

        String utils_docker

        Boolean localize_reads = false

        RuntimeAttr? runtime_attr_localize_bam
        RuntimeAttr? runtime_attr_extract_tags
        RuntimeAttr? runtime_attr_merge_tags
        RuntimeAttr? runtime_attr_extract_contig
        RuntimeAttr? runtime_attr_transfer_tags
        RuntimeAttr? runtime_attr_sort_contig
        RuntimeAttr? runtime_attr_merge_bams
    }

    if (localize_reads) {
        call LocalizeBam {
            input:
                input_bam = aligned_bam,
                input_bai = aligned_bai,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_localize_bam
        }
    }

    File aligned_bam_ = select_first([LocalizeBam.bam, aligned_bam])
    File aligned_bai_ = select_first([LocalizeBam.bai, aligned_bai])

    scatter (path in unaligned_bam_paths) {
        call ExtractMethylationTags {
            input:
                unaligned_bam_path = path,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_tags
        }
    }

    call Helpers.ConcatTsvs as MergeTagsTsvs {
        input:
            tsvs = ExtractMethylationTags.tags_tsv,
            sort_output = false,
            compressed_tsvs = true,
            compressed_output = true,
            prefix = "~{prefix}.all.tags",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_tags
    }

    scatter (contig in contigs) {
        call ExtractContigBam {
            input:
                aligned_bam = aligned_bam_,
                aligned_bai = aligned_bai_,
                contig = contig,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_contig
        }

        call TransferTagsToContig {
            input:
                contig_bam = ExtractContigBam.contig_bam,
                tags_tsv = MergeTagsTsvs.concatenated_tsv,
                prefix = "~{prefix}.~{contig}.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_transfer_tags
        }

        call SortIndexBam {
            input:
                unsorted_bam = TransferTagsToContig.tagged_unsorted_bam,
                prefix = "~{prefix}.~{contig}.sorted",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_sort_contig
        }
    }

    call Helpers.MergeBams {
        input:
            bams = SortIndexBam.tagged_bam,
            bais = SortIndexBam.tagged_bai,
            prefix = "~{prefix}.tagged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_bams
    }

    output {
        File methylation_tagged_bam = MergeBams.merged_bam
        File methylation_tagged_bai = MergeBams.merged_bam_idx
    }
}


task LocalizeBam {
    input {
        File input_bam
        File input_bai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        cp ~{input_bam} localized.bam
        cp ~{input_bai} localized.bam.bai
    >>>

    output {
        File bam = "localized.bam"
        File bai = "localized.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(input_bam, "GB")) + 20,
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

task ExtractMethylationTags {
    input {
        String unaligned_bam_path
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String bam_basename = basename(unaligned_bam_path, ".bam")

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
import os
import pysam

os.environ['HTS_S3_NO_SIGN_REQUEST'] = '1'

with pysam.AlignmentFile("~{unaligned_bam_path}", "rb", check_sq=False) as ubam:
    with gzip.open("~{bam_basename}.tags.tsv.gz", "wt") as out:
        for read in ubam:
            try:
                mm = read.get_tag('MM')
                ml = ','.join(str(v) for v in read.get_tag('ML'))
                out.write(f"{read.query_name}\t{mm}\t{ml}\n")
            except KeyError:
                pass
CODE
    >>>

    output {
        File tags_tsv = "~{bam_basename}.tags.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 25,
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

task ExtractContigBam {
    input {
        File aligned_bam
        File aligned_bai
        String contig
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        aligned_bam: { localization_optional: true }
        aligned_bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        samtools view \
            --threads 3 \
            -h \
            -b \
            -o contig_aligned.bam \
            ~{aligned_bam} \
            ~{contig}
    >>>

    output {
        File contig_bam = "contig_aligned.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 25,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: 4
        memory: 8 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task TransferTagsToContig {
    input {
        File contig_bam
        File tags_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import array
import gzip
import pysam

# First pass: collect the read names present in this contig - only load relevant tags
contig_read_names = set()
with pysam.AlignmentFile("~{contig_bam}", "rb") as abam:
    for read in abam:
        contig_read_names.add(read.query_name)

# Second pass: stream through merged tags, keeping only contig reads
tags_dict = {}
with gzip.open("~{tags_tsv}", 'rt') as f:
    for line in f:
        read_name, mm, ml_str = line.rstrip('\n').split('\t')
        if read_name in contig_read_names:
            tags_dict[read_name] = (mm, [int(v) for v in ml_str.split(',')])

# Third pass: stream through contig BAM, adding tags where relevant
with pysam.AlignmentFile("~{contig_bam}", "rb") as abam:
    with pysam.AlignmentFile("~{prefix}.unsorted.bam", "wb", header=abam.header) as outbam:
        for read in abam:
            if read.query_name in tags_dict:
                mm, ml = tags_dict[read.query_name]
                read.set_tag('MM', mm)
                read.set_tag('ML', array.array('B', ml))
            outbam.write(read)
CODE
    >>>

    output {
        File tagged_unsorted_bam = "~{prefix}.unsorted.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(contig_bam, "GB")) + ceil(size(tags_tsv, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: 1
        memory: 8 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SortIndexBam {
    input {
        File unsorted_bam
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools sort \
            --threads 3 \
            -o "~{prefix}.bam" \
            ~{unsorted_bam}

        samtools index --threads 3 "~{prefix}.bam"
    >>>

    output {
        File tagged_bam = "~{prefix}.bam"
        File tagged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(unsorted_bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: 4
        memory: 8 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
