version 1.0

import "../utils/Structs.wdl"

workflow VcfDist {
    input {
        File vcf_eval
        File? vcf_eval_idx
        File vcf_truth
        File? vcf_truth_idx
        String prefix

        File ref_fa
        File? bed_regions
        String? mode
        Float? threshold
        String? vcfdist_args
        
        String vcfdist_docker

        RuntimeAttr? runtime_attr_vcfdist
    }

    call RunVcfDist {
        input:
            vcf_eval = vcf_eval,
            vcf_eval_idx = vcf_eval_idx,
            vcf_truth = vcf_truth,
            vcf_truth_idx = vcf_truth_idx,
            prefix = prefix,
            ref_fa = ref_fa,
            bed_regions = bed_regions,
            mode = mode,
            threshold = threshold,
            docker = vcfdist_docker,
            runtime_attr_override = runtime_attr_vcfdist,
    }

    output {
        File vcf_dist_summary_vcf = RunVcfDist.summary_vcf
        File vcf_dist_query_tsv = RunVcfDist.query_tsv
        File vcf_dist_truth_tsv = RunVcfDist.truth_tsv
        File vcf_dist_superclusters_tsv = RunVcfDist.superclusters_tsv
        File vcf_dist_precision_recall_summary_tsv = RunVcfDist.precision_recall_summary_tsv
        File vcf_dist_precision_recall_tsv = RunVcfDist.precision_recall_tsv
        File vcf_dist_phase_blocks_tsv = RunVcfDist.phase_blocks_tsv
        File vcf_dist_phasing_summary_tsv = RunVcfDist.phasing_summary_tsv
        File vcf_dist_switchflips_tsv = RunVcfDist.switchflips_tsv
        File vcf_dist_distance_summary = RunVcfDist.distance_summary
        File vcf_dist_distance_tsv = RunVcfDist.distance_tsv
        File vcf_dist_edits_tsv = RunVcfDist.edits_tsv
    }
}

task RunVcfDist {
    input {
        File vcf_eval
        File? vcf_eval_idx
        File vcf_truth
        File? vcf_truth_idx
        String prefix
        File ref_fa
        File? bed_regions
        String? mode
        Float? threshold
        String? vcfdist_args
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if ~{!defined(vcf_eval_idx)}; then
            tabix -p vcf ~{vcf_eval}
        fi

        if ~{!defined(vcf_truth_idx)}; then
            tabix -p vcf ~{vcf_truth}
        fi
        
        mkdir -p results
        
        vcfdist \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{ref_fa} \
            -p "results/~{prefix}." \
            ~{if defined(bed_regions) then "-b " + bed_regions else ""} \
            ~{if defined(mode) then "-m " + mode else ""} \
            ~{if defined(threshold) then "-t " + threshold else ""} \
            ~{if defined(vcfdist_args) then vcfdist_args else ""}
    >>>

    output {
        File summary_vcf = "results/~{prefix}.summary.vcf"
        File query_tsv = "results/~{prefix}.query.tsv"
        File truth_tsv = "results/~{prefix}.truth.tsv"
        File superclusters_tsv = "results/~{prefix}.superclusters.tsv"
        File precision_recall_summary_tsv = "results/~{prefix}.precision-recall-summary.tsv"
        File precision_recall_tsv = "results/~{prefix}.precision-recall.tsv"
        File phase_blocks_tsv = "results/~{prefix}.phase-blocks.tsv"
        File phasing_summary_tsv = "results/~{prefix}.phasing-summary.tsv"
        File switchflips_tsv = "results/~{prefix}.switchflips.tsv"
        File distance_summary = "results/~{prefix}.distance-summary.tsv"
        File distance_tsv = "results/~{prefix}.distance.tsv"
        File edits_tsv = "results/~{prefix}.edits.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(ref_fa, "GB") + size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 10,
        boot_disk_gb: 50,
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