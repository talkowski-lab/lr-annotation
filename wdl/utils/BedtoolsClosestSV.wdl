version 1.0

import "Helpers.wdl"
import "Structs.wdl"

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx
        String prefix

        Int min_sv_length_eval
        Int min_sv_length_truth
        String type_field_eval
        String length_field_eval

        String benchmark_annotations_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_split_eval
        RuntimeAttr? runtime_attr_split_truth
        RuntimeAttr? runtime_attr_compare
        RuntimeAttr? runtime_attr_calculate
        RuntimeAttr? runtime_attr_merge_comparisons
    }

    # Preprocess eval VCF
    call Helpers.SubsetVcfByLength as SubsetEval {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval_idx,
            length_field = length_field_eval,
            min_length = min_sv_length_eval,
            prefix = "~{prefix}.subset_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_eval
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = SubsetEval.subset_vcf,
            vcf_idx = SubsetEval.subset_vcf_idx,
            type_field = type_field_eval,
            length_field = length_field_eval,
            prefix = "~{prefix}.eval.symbolic",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call SplitVcf as SplitEval {
        input:
            vcf = ConvertToSymbolic.processed_vcf,
            vcf_idx = ConvertToSymbolic.processed_vcf_idx,
            split_cpx = false,
            prefix = "~{prefix}.eval",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_split_eval
    }

    # Preprocess truth VCF
    call Helpers.SubsetVcfByLength as SubsetTruth {
        input:
            vcf = vcf_sv_truth,
            vcf_idx = vcf_sv_truth_idx,
            length_field = "SVLEN",
            min_length = min_sv_length_truth,
            prefix = "~{prefix}.subset_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_truth
    }

    call SplitVcf as SplitTruth {
        input:
            vcf = SubsetTruth.subset_vcf,
            vcf_idx = SubsetTruth.subset_vcf_idx,
            split_cpx = true,
            prefix = "~{prefix}.truth",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_split_truth
    }

    # Comparisons by variant type
    call Helpers.BedtoolsClosest as CompareDEL {
        input:
            bed_a = SplitEval.del_bed,
            bed_b = SplitTruth.del_bed,
            prefix = "~{prefix}.DEL",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedSVs as CalcuDEL {
        input:
            input_bed = CompareDEL.output_bed,
            prefix = "~{prefix}.DEL",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareINS {
        input:
            bed_a = SplitEval.ins_bed,
            bed_b = SplitTruth.ins_bed,
            prefix = "~{prefix}.INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuINS {
        input:
            input_bed = CompareINS.output_bed,
            prefix = "~{prefix}.INS",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareDUP {
        input:
            bed_a = SplitEval.dup_bed,
            bed_b = SplitTruth.dup_bed,
            prefix = "~{prefix}.DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedSVs as CalcuDUP {
        input:
            input_bed = CompareDUP.output_bed,
            prefix = "~{prefix}.DUP",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call CollapseRangedToPoint as CollapseTruthDUP {
        input:
            bed = SplitTruth.dup_bed,
            prefix = "~{prefix}.truth.DUP_as_point",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareINS_DUP {
        input:
            bed_a = SplitEval.ins_bed,
            bed_b = CollapseTruthDUP.point_bed,
            prefix = "~{prefix}.INS_DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuINS_DUP {
        input:
            input_bed = CompareINS_DUP.output_bed,
            prefix = "~{prefix}.INS_DUP",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call CollapseRangedToPoint as CollapseEvalDUP {
        input:
            bed = SplitEval.dup_bed,
            prefix = "~{prefix}.eval.DUP_as_point",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareDUP_INS {
        input:
            bed_a = CollapseEvalDUP.point_bed,
            bed_b = SplitTruth.ins_bed,
            prefix = "~{prefix}.DUP_INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuDUP_INS {
        input:
            input_bed = CompareDUP_INS.output_bed,
            prefix = "~{prefix}.DUP_INS",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call PrioritizedConcatComparisons {
        input:
            primary_tsvs = [CalcuDEL.output_comp, CalcuINS.output_comp, CalcuDUP.output_comp],
            secondary_tsvs = [CalcuINS_DUP.output_comp, CalcuDUP_INS.output_comp],
            prefix = "~{prefix}.comparison",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    call CreateBedtoolsAnnotationTsv {
        input:
            truvari_unmatched_vcf = SubsetEval.subset_vcf,
            truvari_unmatched_vcf_idx = SubsetEval.subset_vcf_idx,
            closest_bed = PrioritizedConcatComparisons.merged_tsv,
            prefix = "~{prefix}.bedtools_closest_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    output {
        File annotation_tsv = CreateBedtoolsAnnotationTsv.annotation_tsv
    }
}

task SplitVcf {
    input {
        File vcf
        File vcf_idx
        Boolean split_cpx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed \
            -i SVTYPE \
            -i SVLEN \
            ~{if split_cpx then "--split-cpx" else ""} \
            ~{vcf} \
            tmp.bed

        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed

        set +o pipefail

        head -1 ~{prefix}.bed > header

        set -o pipefail

        cat header <(awk '$5 == "DEL"' ~{prefix}.bed) > ~{prefix}.DEL.bed
        cat header <(awk '$5 == "DUP"' ~{prefix}.bed) > ~{prefix}.DUP.bed
        cat header <(awk '$5 ~ /^INS/' ~{prefix}.bed) > ~{prefix}.INS.bed
        cat header <(awk '$5 == "INV" || $5 == "CPX"' ~{prefix}.bed) > ~{prefix}.INV_CPX.bed
        cat header <(awk '$5 == "BND" || $5 == "CTX"' ~{prefix}.bed) > ~{prefix}.BND_CTX.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del_bed = "~{prefix}.DEL.bed"
        File dup_bed = "~{prefix}.DUP.bed"
        File ins_bed = "~{prefix}.INS.bed"
        File inv_bed = "~{prefix}.INV_CPX.bed"
        File bnd_bed = "~{prefix}.BND_CTX.bed"
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

task SelectMatchedSVs {
    input {
        File input_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        Rscript /opt/gnomad-lr/scripts/benchmark/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task SelectMatchedINSs {
    input {
        File input_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        Rscript /opt/gnomad-lr/scripts/benchmark/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CollapseRangedToPoint {
    input {
        File bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        head -1 ~{bed} > ~{prefix}.bed

        awk 'NR>1 {OFS="\t"; $3=$2+1; print}' ~{bed} >> ~{prefix}.bed
    >>>

    output {
        File point_bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed, "GB")) + 5,
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

task PrioritizedConcatComparisons {
    input {
        Array[File] primary_tsvs
        Array[File] secondary_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " primary_tsvs} > ~{prefix}.comparison

        grep -v "query_svid" ~{prefix}.comparison | cut -f1 | sort -u > seen_ids.txt

        for f in ~{sep=" " secondary_tsvs}; do
            awk -F'\t' 'NR==FNR{seen[$1]=1; next} !($1 in seen)' seen_ids.txt "$f" >> ~{prefix}.comparison
        done
    >>>

    output {
        File merged_tsv = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CreateBedtoolsAnnotationTsv {
    input {
        File truvari_unmatched_vcf
        File truvari_unmatched_vcf_idx
        File closest_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export LC_ALL=C

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{truvari_unmatched_vcf} \
        | sort -k5,5 > vcf_info.tsv

        grep -v "query_svid" ~{closest_bed} \
            | awk -F'\t' '$1 != "" {print $1"\t"$2}' \
            | sort -k1,1 > matched_ids.tsv

        join \
            -1 5 \
            -2 1 \
            -t $'\t' \
            vcf_info.tsv \
            matched_ids.tsv \
            > joined.tsv

        awk -F'\t' 'BEGIN{OFS="\t"} {
            print $2, $3, $4, $5, $1, "BEDTOOLS_CLOSEST", $6, "SV"
        }' joined.tsv \
        | sort -k1,1V -k2,2n > ~{prefix}.bedtools_matched.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.bedtools_matched.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(truvari_unmatched_vcf, "GB") + size(closest_bed, "GB")) + 5,
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
