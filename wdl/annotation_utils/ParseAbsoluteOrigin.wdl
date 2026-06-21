version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ParseAbsoluteOrigin {
    input {
        File vcf
        File vcf_idx
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_parse
        RuntimeAttr? runtime_attr_concat_tsvs
        RuntimeAttr? runtime_attr_annotate
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = "~{prefix}.sharded",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [vcf]])
    Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx]])

    scatter (i in range(length(vcfs_to_process))) {
        call Helpers.SubsetVcfByArgs as SubsetFlankOrigin {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                include_args = 'INFO/ORIGIN~"^flank_"',
                prefix = "~{prefix}.shard_~{i}.flank_origin_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call ParseShardOrigin {
            input:
                vcf = SubsetFlankOrigin.subset_vcf,
                vcf_idx = SubsetFlankOrigin.subset_vcf_idx,
                prefix = "~{prefix}.shard_~{i}.origin_updates",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_parse
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = ParseShardOrigin.origin_tsv,
            sort_output = true,
            prefix = "~{prefix}.origin_updates",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    call AnnotateAbsoluteOrigin {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            origin_tsv = ConcatTsvs.concatenated_tsv,
            prefix = "~{prefix}.absolute_origin",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_annotate
    }

    output {
        File absolute_origin_vcf = AnnotateAbsoluteOrigin.annotated_vcf
        File absolute_origin_vcf_idx = AnnotateAbsoluteOrigin.annotated_vcf_idx
    }
}

task ParseShardOrigin {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

vcf_in = pysam.VariantFile("~{vcf}")
with open("~{prefix}.tsv", "w") as out:
    for record in vcf_in:
        origin = record.info.get("ORIGIN")
        if origin is None:
            continue
        origin_values = origin if isinstance(origin, (list, tuple)) else (origin,)
        flank_origin = next((v for v in origin_values if v and v.startswith("flank_")), None)
        if flank_origin is None:
            continue
        body, strand = flank_origin.rsplit("_", 1)
        _, coord_part = body.rsplit("_", 1)
        abs_chrom, _, local_range = coord_part.split(":")
        local_start_str, local_end_str = local_range.split("-")
        local_start = int(local_start_str)
        local_end = int(local_end_str)
        allele_length = record.info["allele_length"]
        if isinstance(allele_length, (list, tuple)):
            allele_length = allele_length[0]
        alt_len = abs(int(allele_length))
        flank_start = max(1, record.pos - alt_len - 100)
        abs_start = flank_start + local_start
        abs_end = flank_start + local_end
        new_origin = f"{abs_chrom}:{abs_start}-{abs_end}_{strand}"
        alts = ",".join(record.alts) if record.alts else "."
        rid = record.id if record.id else "."
        out.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{alts}\t{rid}\t{new_origin}\n")
vcf_in.close()
CODE
    >>>

    output {
        File origin_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task AnnotateAbsoluteOrigin {
    input {
        File vcf
        File vcf_idx
        File origin_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bgzip -c ~{origin_tsv} > origin_updates.tsv.gz
        tabix -s1 -b2 -e2 origin_updates.tsv.gz

        bcftools annotate \
            -a origin_updates.tsv.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/ORIGIN \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(origin_tsv, "GB")) + 10,
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
