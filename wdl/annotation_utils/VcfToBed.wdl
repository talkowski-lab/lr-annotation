version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow VcfToBed {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Boolean include_samples = true
        Boolean include_filters = true
        Boolean split_bnd = false
        Boolean split_cpx = false

        String gatk_sv_lr_docker
        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_concat
    }

    scatter (i in range(length(contigs))) {
        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contigs[i]}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [vcfs[i]]])

        scatter (j in range(length(shard_vcfs))) {
            call SvtkVcfToBed {
                input:
                    vcf = shard_vcfs[j],
                    include_samples = include_samples,
                    include_filters = include_filters,
                    split_bnd = split_bnd,
                    split_cpx = split_cpx,
                    prefix = "~{prefix}.~{contigs[i]}.shard_~{j}",
                    docker = gatk_sv_lr_docker,
                    runtime_attr_override = runtime_attr_vcf2bed
            }
        }
    }

    call Helpers.ConcatTsvs as ConcatBedShards {
        input:
            tsvs = flatten(SvtkVcfToBed.bed),
            sort_output = false,
            preserve_header = true,
            compressed_output = true,
            prefix = "~{prefix}.vcf2bed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File bed = ConcatBedShards.concatenated_tsv
    }
}

task SvtkVcfToBed {
    input {
        File vcf
        Boolean include_samples
        Boolean include_filters
        Boolean split_bnd
        Boolean split_cpx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed \
            --info ALL \
            ~{if include_samples then "" else "--no-samples"} \
            ~{if include_filters then "--include-filters" else ""} \
            ~{if split_bnd then "--split-bnd" else ""} \
            ~{if split_cpx then "--split-cpx" else ""} \
            ~{vcf} \
            ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
