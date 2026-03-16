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

        RuntimeAttr? runtime_attr_extract_tags
        RuntimeAttr? runtime_attr_transfer_tags
        RuntimeAttr? runtime_attr_merge_bams
    }

    scatter (path in unaligned_bam_paths) {
        call ExtractMethylationTags {
            input:
                unaligned_bam_path = path,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_tags
        }
    }

    scatter (contig in contigs) {
        call TransferTagsToContig {
            input:
                aligned_bam = aligned_bam,
                aligned_bai = aligned_bai,
                contig = contig,
                tags_tsvs = ExtractMethylationTags.tags_tsv,
                prefix = "~{prefix}.~{contig}.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_transfer_tags
        }
    }

    call Helpers.MergeBams {
        input:
            bams = TransferTagsToContig.tagged_bam,
            bais = TransferTagsToContig.tagged_bai,
            prefix = "~{prefix}.tagged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_bams
    }

    output {
        File methylation_tagged_bam = MergeBams.merged_bam
        File methylation_tagged_bai = MergeBams.merged_bam_idx
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
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 30,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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

task TransferTagsToContig {
    input {
        File aligned_bam
        File aligned_bai
        String contig
        Array[File] tags_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "Subsetting aligned BAM to ~{contig}..."
        samtools view \
            -h -b \
            -@ 3 \
            ~{aligned_bam} \
            ~{contig} \
            -o contig_aligned.bam

        echo "Transferring MM and ML tags to ~{contig}..."
        python3 <<CODE
import array
import gzip
import pysam

with open("~{write_lines(tags_tsvs)}") as f:
    tsv_paths = f.read().strip().split('\n')

tags_dict = {}
for tsv_path in tsv_paths:
    with gzip.open(tsv_path, 'rt') as f:
        for line in f:
            read_name, mm, ml_str = line.rstrip('\n').split('\t')
            tags_dict[read_name] = (mm, [int(v) for v in ml_str.split(',')])

with pysam.AlignmentFile("contig_aligned.bam", "rb") as abam:
    with pysam.AlignmentFile("~{prefix}.unsorted.bam", "wb", header=abam.header) as outbam:
        for read in abam:
            if read.query_name in tags_dict:
                mm, ml = tags_dict[read.query_name]
                read.set_tag('MM', mm)
                read.set_tag('ML', array.array('B', ml))
            outbam.write(read)
CODE

        echo "Sorting and indexing..."
        samtools sort -@ 3 -o "~{prefix}.bam" "~{prefix}.unsorted.bam"
        samtools index -@ 3 "~{prefix}.bam"
    >>>

    output {
        File tagged_bam = "~{prefix}.bam"
        File tagged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 64,
        disk_gb: ceil(size(aligned_bam, "GB") * 0.15) * 2 + ceil(size(tags_tsvs, "GB")) + 30,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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
