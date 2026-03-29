version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CallInsRemap {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Array[File] ref_bwa_indices
        String prefix
        Int? records_per_shard
        Int? minlength
        Int? maxlength
        Int? mm2_threshold
        Float? cov_threshold

        String utils_docker = "us.gcr.io/broad-dsp-lrma/lr-basic:latest"

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
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
        call InsRemap {
            input:
                vcf = vcfs_to_process[shard_idx],
                vcf_idx = vcf_idxs_to_process[shard_idx],
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                ref_bwa_indices = ref_bwa_indices,
                prefix = "~{prefix}.shard_~{shard_idx}",
                minlength = minlength,
                maxlength = maxlength,
                mm2_threshold = mm2_threshold,
                cov_threshold = cov_threshold
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatTsvs as ConcatShards {
            input:
                tsvs = InsRemap.remap_tsv,
                prefix = "~{prefix}.remap_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    File final_remap_tsv = select_first([ConcatShards.concatenated_tsv, InsRemap.remap_tsv[0]])
    
    output {
        File remap_tsv = final_remap_tsv
    }
}

task InsRemap {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Array[File] ref_bwa_indices
        String prefix
        Int minlength=50
        Int maxlength=10000000 #don't attempt remap for INS over 10 Mb
        Int mm2_threshold=5000 #at this length, switch from bwa to minimap2 for alignment
        Float cov_threshold=0.8
        
        RuntimeAttr? runtime_attr_override
    }
   
    Int disk_size = ceil(size(vcf, "GB") + size(ref_fa, "GB")) * 3 + 50

    command <<<
        set -euxo pipefail
        echo "starting"
        # how many cpus are available?
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "moving refs"
        mkdir ref_files
        mv ~{ref_fa} ref_files/
        mv ~{ref_fai} ref_files/
        
        for x in ~{sep=' ' ref_bwa_indices}; do
            mv $x ref_files/
        done

        ref_fa_basename=$(basename ~{ref_fa})
        echo "running truvari"
        truvari anno remap -r ref_files/$ref_fa_basename -o ~{prefix}.remapped.vcf.gz \
            --min-length ~{minlength} --max-length ~{maxlength} --mm2-threshold ~{mm2_threshold} \
            --threads ${N_THREADS} --cov-threshold ~{cov_threshold} ~{vcf}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%remap_classification\t%remap_coords\t%remap_ori\t%remap_perc\t%remap_segments\n' \
            ~{prefix}.remapped.vcf.gz  > ~{prefix}.tsv

        tabix -f -p vcf ~{prefix}.remapped.vcf.gz
    >>>
    
    output {
        File remapped_vcf = "~{prefix}.remapped.vcf.gz"
        File remapped_vcf_idx = "~{prefix}.remapped.vcf.gz.tbi"
        File remap_tsv = "~{prefix}.tsv"
    }

     #########################
     RuntimeAttr default_attr = object {
         cpu_cores:          8,
         mem_gb:             16,
         disk_gb:            disk_size,
         boot_disk_gb:       10,
         preemptible_tries:  2,
         max_retries:        0,
         docker:             "quay.io/ymostovoy/lr-remap:latest"
     }
     RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
     runtime {
         cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
         memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
         disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
         bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
         preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
         maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
         docker:                 select_first([runtime_attr.docker,            default_attr.docker])
     }
}