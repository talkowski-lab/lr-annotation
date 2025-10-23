version 1.0

import "../general/Structs.wdl"

workflow ShardAndComputeBenchmarks {
    input {
        File final_vcf
        File final_vcf_index
        File matched_ids_tsv
        File truth_tsv_snv
        File truth_tsv_sv
        File truth_vep_header
        String contig
        String prefix
        String pipeline_docker
        Int variants_per_shard
        String? skip_vep_categories

        RuntimeAttr? runtime_attr_shard_matched_eval
        RuntimeAttr? runtime_attr_compute_shard_benchmarks
        RuntimeAttr? runtime_attr_merge_shard_benchmarks
    }

    call ShardMatchedEval {
        input:
            matched_ids_tsv = matched_ids_tsv,
            variants_per_shard = variants_per_shard,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_shard_matched_eval
    }

    scatter (shard_idx in range(length(ShardMatchedEval.shard_id_lists))) {
        call ComputeShardBenchmarks {
            input:
                final_vcf = final_vcf,
                final_vcf_index = final_vcf_index,
                matched_shard_tsv = ShardMatchedEval.shard_id_lists[shard_idx],
                truth_tsv_snv = truth_tsv_snv,
                truth_tsv_sv = truth_tsv_sv,
                truth_vep_header = truth_vep_header,
                contig = contig,
                shard_label = "~{shard_idx}",
                prefix = prefix,
                skip_vep_categories = skip_vep_categories,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_compute_shard_benchmarks
        }
    }

    call MergeShardBenchmarks {
        input:
            af_pair_tsvs = select_all(ComputeShardBenchmarks.af_pairs_tsv),
            vep_pair_tsvs = select_all(ComputeShardBenchmarks.vep_pairs_tsv),
            truth_vep_header = truth_vep_header,
            contig = contig,
            prefix = prefix,
            skip_vep_categories = skip_vep_categories,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge_shard_benchmarks
    }

    output {
        File plot_tarball = MergeShardBenchmarks.plot_tarball
    }
}

task ShardMatchedEval {
    input {
        File matched_ids_tsv
        Int variants_per_shard
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        mkdir -p shards
        zcat ~{matched_ids_tsv} | awk 'BEGIN{c=0;f=0} {print > sprintf("shards/ids.%06d.tsv", int(c/~{variants_per_shard})) ; c++} END{ }'
        ls shards/ids.*.tsv | while read f; do bgzip -f "$f"; done
    >>>

    output {
        Array[File] shard_id_lists = glob("shards/*.tsv.gz")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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

task ComputeShardBenchmarks {
    input {
        File final_vcf
        File final_vcf_index
        File matched_shard_tsv
        File truth_tsv_snv
        File truth_tsv_sv
        File truth_vep_header
        String contig
        String shard_label
        String prefix
        String? skip_vep_categories
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_benchmarks_shard.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --final_vcf ~{final_vcf} \
            --matched_shard_tsv ~{matched_shard_tsv} \
            --truth_tsv_snv ~{truth_tsv_snv} \
            --truth_tsv_sv ~{truth_tsv_sv} \
            --truth_vep_header ~{truth_vep_header} \
            --shard_label ~{shard_label} \
            ~{"--skip_vep_categories " + skip_vep_categories}
    >>>

    output {
        File af_pairs_tsv = "~{prefix}.~{contig}.shard_~{shard_label}.af_pairs.tsv.gz"
        File vep_pairs_tsv = "~{prefix}.~{contig}.shard_~{shard_label}.vep_pairs.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(final_vcf, "GB")) * 2 + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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

task MergeShardBenchmarks {
    input {
        Array[File] af_pair_tsvs
        Array[File] vep_pair_tsvs
        File truth_vep_header
        String contig
        String prefix
        String? skip_vep_categories
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/benchmark/merge_benchmarks_from_pairs.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --af_pair_tsvs ~{sep=',' af_pair_tsvs} \
            --vep_pair_tsvs ~{sep=',' vep_pair_tsvs} \
            --truth_vep_header ~{truth_vep_header} \
            ~{"--skip_vep_categories " + skip_vep_categories}
    >>>

    output {
        File plot_tarball = "~{prefix}.benchmarks.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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
