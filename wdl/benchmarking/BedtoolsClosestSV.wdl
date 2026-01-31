version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx
        String prefix
        
        String benchmark_annotations_docker
        String utils_docker
        
        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_split_eval
        RuntimeAttr? runtime_attr_split_truth
        RuntimeAttr? runtime_attr_compare
        RuntimeAttr? runtime_attr_calculate
        RuntimeAttr? runtime_attr_merge_comparisons
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval_idx,
            prefix = "~{prefix}.eval.symbolic",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call SplitQueryVcf as SplitEval {
        input:
            vcf = ConvertToSymbolic.processed_vcf,
            vcf_idx = ConvertToSymbolic.processed_vcf_idx,
            type_field = "allele_type",
            length_field = "allele_length",
            prefix = "~{prefix}.eval",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_split_eval
    }

    call SplitQueryVcf as SplitTruth {
        input:
            vcf = vcf_sv_truth,
            vcf_idx = vcf_sv_truth_idx,
            type_field = "allele_type",
            length_field = "SVLEN",
            prefix = "~{prefix}.truth",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_split_truth
    }

    call Helpers.BedtoolsClosest as CompareDEL {
        input:
            bed_a = SplitEval.del_bed,
            bed_b = SplitTruth.del_bed,
            allele_type = "DEL",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }
    
    call SelectMatchedSVs as CalcuDEL {
        input:
            input_bed = CompareDEL.output_bed,
            allele_type = "DEL",
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareINS {
        input:
            bed_a = SplitEval.ins_bed,
            bed_b = SplitTruth.ins_bed,
            allele_type = "INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }
    
    call SelectMatchedINSs as CalcuINS {
        input:
            input_bed = CompareINS.output_bed,
            docker = benchmark_annotations_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = [CalcuDEL.output_comp, CalcuINS.output_comp],
            prefix = "~{prefix}.comparison",
            skip_sort = true,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    call CreateBedtoolsAnnotationTsv {
        input:
            truvari_unmatched_vcf = vcf_eval,
            truvari_unmatched_vcf_idx = vcf_eval_idx,
            closest_bed = ConcatTsvs.concatenated_tsv,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    output {
        File closest_bed = ConcatTsvs.concatenated_tsv
        File annotation_tsv = CreateBedtoolsAnnotationTsv.annotation_tsv
    }
}

task SplitQueryVcf {
    input {
        File vcf
        File vcf_idx
        String type_field
        String length_field
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed \
            -i ~{type_field} \
            -i ~{length_field} \
            ~{vcf} \
            tmp.bed
        
        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed

        set +o pipefail

        head -1 ~{prefix}.bed > header

        set -o pipefail
        
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed )> ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed )> ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' ~{prefix}.bed )> ~{prefix}.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' ~{prefix}.bed )> ~{prefix}.INV_CPX.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' ~{prefix}.bed )> ~{prefix}.BND_CTX.bed
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
        String allele_type
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_bed, ".bed")

    command <<<
        set -euo pipefail

        Rscript /opt/gnomad-lr/scripts/benchmark/R_scripts/R1.bedtools_closest_CNV.R \
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
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_bed, ".bed")

    command <<<
        set -euo pipefail
        
        Rscript /opt/gnomad-lr/scripts/benchmark/R_scripts/R2.bedtools_closest_INS.R \
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
        }' joined.tsv > ~{prefix}.bedtools_matched.tsv
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
