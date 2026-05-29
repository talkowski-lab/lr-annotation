version 1.0

workflow AnnotateVcf {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        String prefix

        Int? records_per_shard

        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(records_per_shard)) {
        call ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [vcf]])
    Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx]])

    scatter (shard_idx in range(length(vcfs_to_process))) {
        call AnnotateSequentially as AnnotateVcfShard {
            input:
                vcf = vcfs_to_process[shard_idx],
                vcf_idx = vcf_idxs_to_process[shard_idx],
                annotations_tsvs = annotations_tsvs,
                info_names = info_names,
                info_descriptions = info_descriptions,
                info_types = info_types,
                info_numbers = info_numbers,
                prefix = "~{prefix}.annotated.shard_~{shard_idx}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    if (defined(records_per_shard)) {
        call ConcatBcfs {
            input:
                bcfs = AnnotateVcfShard.annotated_vcf,
                bcf_idxs = AnnotateVcfShard.annotated_vcf_idx,
                prefix = prefix + ".annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotated_vcf = select_first([ConcatBcfs.concat_bcf, AnnotateVcfShard.annotated_vcf[0]])
        File annotated_vcf_idx = select_first([ConcatBcfs.concat_bcf_idx, AnnotateVcfShard.annotated_vcf_idx[0]])
    }
}

task ShardVcfByRecords {
    input {
        File vcf
        File vcf_idx
        Int records_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(10 + 3 * size(vcf, "GB")),
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        mkdir -p shards

        bcftools +scatter ~{vcf} \
            -o shards \
            -O z \
            --prefix ~{prefix}. \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -n ~{records_per_shard}

        for shard in shards/*.vcf.gz; do
            tabix -p vcf "$shard"
        done
    >>>

    output {
        Array[File] shards = glob("shards/*.vcf.gz")
        Array[File] shard_idxs = glob("shards/*.vcf.gz.tbi")
    }

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

task AnnotateSequentially {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
import sys

info_names = [line.strip().split('\t') for line in open('~{write_tsv(info_names)}')]
info_descriptions = [line.strip().split('\t') for line in open('~{write_tsv(info_descriptions)}')]
info_types = [line.strip().split('\t') for line in open('~{write_tsv(info_types)}')]
info_numbers = [line.strip().split('\t') for line in open('~{write_tsv(info_numbers)}')]

if len(info_names) != len(info_descriptions) or len(info_names) != len(info_types) or len(info_names) != len(info_numbers):
    sys.stderr.write("Error: All info arrays must have the same length\n")
    sys.exit(1)

for i, (names, descs, types, numbers) in enumerate(zip(info_names, info_descriptions, info_types, info_numbers)):
    if len(names) != len(descs) or len(names) != len(types) or len(names) != len(numbers):
        sys.stderr.write(f"Error: info arrays at index {i} must have the same length\n")
        sys.exit(1)
    
    with open(f"header_{i}.txt", "w") as f:
        for name, desc, type_val, number in zip(names, descs, types, numbers):
            f.write(f'##INFO=<ID={name},Number={number},Type={type_val},Description="{desc}">\n')
    
    column_spec = ','.join(['CHROM', 'POS', 'REF', 'ALT', '~ID'] + [f'INFO/{name}' for name in names])
    with open(f"columns_{i}.txt", "w") as f:
        f.write(column_spec)

with open("num_tsvs.txt", "w") as f:
    f.write(str(len(info_names)))
EOF

        current_vcf="~{vcf}"
        
        i=0
        for tsv_file in ~{sep=' ' annotations_tsvs}; do
            sort -k1,1 -k2,2n "$tsv_file" |bgzip -c > "annotations_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"
            
            COLUMN_SPEC=$(cat "columns_${i}.txt")
            
            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -h "header_${i}.txt" \
                -c "$COLUMN_SPEC" \
                -Oz -o "temp_${i}.vcf.gz" \
                "$current_vcf"
            current_vcf="temp_${i}.vcf.gz"
            
            i=$((i + 1))
        done
        
        # hardcoding this line just for MEIs!
        bcftools view -Ob -i 'MEI_type!="."' "$current_vcf" -o ~{prefix}.MEIs.bcf
        bcftools index -f ~{prefix}.MEIs.bcf
    >>>

    output {
        File annotated_vcf = "~{prefix}.MEIs.bcf"
        File annotated_vcf_idx = "~{prefix}.MEIs.bcf.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(annotations_tsvs, "GB")) + 20,
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

task ConcatBcfs {
    input {
        Array[File] bcfs
        Array[File] bcf_idxs
        String prefix = "concat"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bcfs, "GB")) + 5,
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        BCFS_FILE="~{write_lines(bcfs)}"

        bcftools concat \
            --allow-overlaps \
            --file-list ${BCFS_FILE} \
            -Ob -o "~{prefix}.MEIs.bcf"

        bcftools index -f "~{prefix}.MEIs.bcf"
    >>>

    output {
        File concat_bcf = "~{prefix}.MEIs.bcf"
        File concat_bcf_idx = "~{prefix}.MEIs.bcf.csi"
    }

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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}