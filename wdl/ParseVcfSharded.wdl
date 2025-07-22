version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow ParseVcfSharded {
    input {
        File vcf
        File vcf_index
        String prefix
        Int variants_per_shard
        String pipeline_docker

        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_parse_vcf
        RuntimeAttr? runtime_attr_concat_tsv
    }

    call Helpers.SplitVcfIntoShards {
        input:
            input_vcf = vcf,
            input_vcf_index = vcf_index,
            variants_per_shard = variants_per_shard,
            output_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    scatter (shard in zip(SplitVcfIntoShards.split_vcfs, SplitVcfIntoShards.split_vcf_indexes)) {
        String shard_prefix = basename(shard.left, ".vcf.gz")
        call ParseVcfToTsv {
            input:
                vcf = shard.left,
                prefix = shard_prefix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_parse_vcf
        }
    }

    call Helpers.ConcatFiles as ConcatTsvs {
        input:
            files = ParseVcfToTsv.output_tsv,
            outfile_name = "~{prefix}.tsv",
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_concat_tsv
    }

    output {
        File output_tsv = ConcatTsvs.concatenated_file
    }
}

task ParseVcfToTsv {
    input {
        File vcf
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/parse_vcf_to_tsv.py ~{vcf} ~{prefix}.tsv
    >>>

    output {
        File output_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB") * 2) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
} 