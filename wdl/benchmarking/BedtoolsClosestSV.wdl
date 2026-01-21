version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_sv_truth
        File vcf_sv_truth_index
        String prefix
        String bedtools_closest_docker
        
        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_split_eval
        RuntimeAttr? runtime_attr_split_truth
        RuntimeAttr? runtime_attr_compare_del
        RuntimeAttr? runtime_attr_calcu_del
        RuntimeAttr? runtime_attr_compare_dup
        RuntimeAttr? runtime_attr_calcu_dup
        RuntimeAttr? runtime_attr_compare_ins
        RuntimeAttr? runtime_attr_calcu_ins
        RuntimeAttr? runtime_attr_compare_inv
        RuntimeAttr? runtime_attr_calcu_inv
        RuntimeAttr? runtime_attr_compare_bnd
        RuntimeAttr? runtime_attr_calcu_bnd
        RuntimeAttr? runtime_attr_merge_comparisons
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval_index,
            prefix = "~{prefix}.eval.symbolic",
            drop_genotypes = true,
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call Helpers.SplitQueryVcf as SplitEval {
        input:
            vcf = ConvertToSymbolic.processed_vcf,
            prefix = "~{prefix}.eval",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_split_eval
    }

    call Helpers.SplitQueryVcf as SplitTruth {
        input:
            vcf = vcf_sv_truth,
            vcf_idx = vcf_sv_truth_index,
            prefix = "~{prefix}.truth",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_split_truth
    }

    call Helpers.BedtoolsClosest as CompareDEL {
        input:
            bed_a = SplitEval.del_bed,
            bed_b = SplitTruth.del_bed,
            svtype = "DEL",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_compare_del
    }
    
    call SelectMatchedSVs as CalcuDEL {
        input:
            input_bed = CompareDEL.output_bed,
            svtype = "DEL",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_calcu_del
    }

    call Helpers.BedtoolsClosest as CompareDUP {
        input:
            bed_a = SplitEval.dup_bed,
            bed_b = SplitTruth.dup_bed,
            svtype = "DUP",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_compare_dup
    }
    
    call SelectMatchedSVs as CalcuDUP {
        input:
            input_bed = CompareDUP.output_bed,
            svtype = "DUP",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_calcu_dup
    }

    call Helpers.BedtoolsClosest as CompareINS {
        input:
            bed_a = SplitEval.ins_bed,
            bed_b = SplitTruth.ins_bed,
            svtype = "INS",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_compare_ins
    }
    
    call SelectMatchedINSs as CalcuINS {
        input:
            input_bed = CompareINS.output_bed,
            svtype = "INS",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_calcu_ins
    }

    call Helpers.BedtoolsClosest as CompareINV {
        input:
            bed_a = SplitEval.inv_bed,
            bed_b = SplitTruth.inv_bed,
            svtype = "INV",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_compare_inv
    }
    
    call SelectMatchedSVs as CalcuINV {
        input:
            input_bed = CompareINV.output_bed,
            svtype = "INV",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_calcu_inv
    }

    call Helpers.BedtoolsClosest as CompareBND {
        input:
            bed_a = SplitEval.bnd_bed,
            bed_b = SplitTruth.bnd_bed,
            svtype = "BND",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_compare_bnd
    }
    
    call SelectMatchedINSs as CalcuBND {
        input:
            input_bed = CompareBND.output_bed,
            svtype = "BND",
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_calcu_bnd
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = [CalcuDEL.output_comp, CalcuDUP.output_comp, CalcuINS.output_comp, CalcuINV.output_comp, CalcuBND.output_comp],
            prefix = "~{prefix}.comparison",
            skip_sort = true,
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    call CreateBedtoolsAnnotationTsv {
        input:
            truvari_unmatched_vcf = vcf_eval,
            truvari_unmatched_vcf_index = vcf_eval_index,
            closest_bed = ConcatTsvs.concatenated_tsv,
            prefix = prefix,
            docker = bedtools_closest_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    output {
        File closest_bed = ConcatTsvs.concatenated_tsv
        File annotation_tsv = CreateBedtoolsAnnotationTsv.annotation_tsv
    }
} 

task SelectMatchedSVs {
    input {
        File input_bed
        String svtype
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

task SelectMatchedINSs {
    input {
        File input_bed
        String svtype
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_bed, ".bed")

    command <<<
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

task CreateBedtoolsAnnotationTsv {
    input {
        File truvari_unmatched_vcf
        File truvari_unmatched_vcf_index
        File closest_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -F'\t' 'NR>1 && $NF!="NA" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tBEDTOOLS_CLOSEST\t"$NF"\tSV"}' \
            ~{closest_bed} > ~{prefix}.bedtools_matched.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.bedtools_matched.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
