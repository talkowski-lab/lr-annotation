version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow TruvariRemap {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        File ref_fa
        Array[File] ref_bwa_indices

        Int min_length
        Int max_length
        Int mm2_threshold
        Float cov_threshold

        String utils_docker

        RuntimeAttr? runtime_attr_ins_remap
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call InsRemap {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                ref_fa = ref_fa,
                ref_bwa_indices = ref_bwa_indices,
                min_length = min_length,
                max_length = max_length,
                mm2_threshold = mm2_threshold,
                cov_threshold = cov_threshold,
                prefix = "~{prefix}.~{contig}.remapped",
                runtime_attr_override = runtime_attr_ins_remap
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = InsRemap.remap_annotations_tsv,
            prefix = "~{prefix}.remap_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File annotations_tsv_remap = ConcatTsvs.concatenated_tsv
    }
}

task InsRemap {
    input {
        File vcf
        File vcf_idx
        String contig
        File ref_fa
        Array[File] ref_bwa_indices
        Int min_length
        Int max_length
        Int mm2_threshold
        Float cov_threshold
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        bcftools view \
            -r ~{contig} \
            -Oz -o ~{contig}.vcf.gz \
            ~{vcf}
        
        tabix ~{contig}.vcf.gz

        mkdir ref_files
        mv ~{ref_fa} ref_files/
        for x in ~{sep=' ' ref_bwa_indices}; do
            mv $x ref_files/
        done

        ref_fa_basename=$(basename ~{ref_fa})

        truvari anno remap \
            -r ref_files/$ref_fa_basename \
            -o ~{prefix}.vcf.gz \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --mm2-threshold ~{mm2_threshold} \
            --threads ${N_THREADS} \
            --cov-threshold ~{cov_threshold} \
            ~{contig}.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/remap_classification\t%INFO/remap_coords\t%INFO/remap_ori\t%INFO/remap_perc\n' \
            ~{prefix}.vcf.gz \
        | awk -F'\t' '($6 != "." || $7 != "." || $8 != "." || $9 != ".")' \
        > ~{prefix}.tsv
    >>>
    
    output {
        File remap_annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fa, "GB")) * 3 + 50,
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
        docker: "quay.io/ymostovoy/lr-remap:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
