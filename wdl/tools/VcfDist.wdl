version 1.0

import "../utils/Structs.wdl"

workflow RunVcfdistBenchmark {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix

        File ref_fa
        File? bed_regions
        String? mode
        Float? threshold
        
        String vcfdist_docker
    }

    call Vcfdist {
        input:
            vcf_eval = vcf_eval,
            vcf_eval_idx = vcf_eval_idx,
            vcf_truth = vcf_truth,
            vcf_truth_idx = vcf_truth_idx,
            ref_fasta = ref_fasta,
            bed_regions = bed_regions,
            prefix = sample_prefix,
            docker = vcfdist_docker
    }

    output {
        File vcfdist_pr_summary = Vcfdist.precision_recall_summary
        File vcfdist_pr_tsv = Vcfdist.precision_recall_tsv
        File vcfdist_dist_summary = Vcfdist.distance_summary
        File vcfdist_dist_tsv = Vcfdist.distance_tsv
    }
}

task task RunVcfDist {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix
        File ref_fa
        File? bed_regions
        String? mode
        Float? threshold
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        mkdir -p vcfdist_out
        
        vcfdist \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{ref_fa} \
            -o "vcfdist_out" \
            -p ~{prefix} \
            ~{if defined(bed_regions) then "-b " + bed_regions else ""} \
            ~{if defined(mode) then "-m " + mode else ""} \
            ~{if defined(threshold) then "-t " + threshold else ""}
    <<<

    output {
        File precision_recall_summary = "vcfdist_out/~{prefix}.precision-recall-summary.tsv"
        File precision_recall_tsv = "vcfdist_out/~{prefix}.precision-recall.tsv"
        File distance_summary = "vcfdist_out/~{prefix}.distance-summary.tsv"
        File distance_tsv = "vcfdist_out/~{prefix}.distance.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth)) + 10,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}