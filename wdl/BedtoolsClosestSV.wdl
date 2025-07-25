version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_truth
        String prefix
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_select_matched_svs
        RuntimeAttr? runtime_attr_concat
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval + ".tbi",
            prefix = "~{prefix}.eval.symbolic",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call Helpers.SplitQueryVcf as SplitEval {
        input:
            vcf = ConvertToSymbolic.processed_vcf,
            prefix = "~{prefix}.eval",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    call Helpers.SplitQueryVcf as SplitTruth {
        input:
            vcf = vcf_truth,
            prefix = "~{prefix}.truth",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    call Helpers.BedtoolsClosest as CompareDEL {
        input:
            bed_a = SplitEval.del_bed,
            bed_b = SplitTruth.del_bed,
            svtype = "DEL",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }
    
    call Helpers.SelectMatchedSVs as CalcuDEL {
        input:
            input_bed = CompareDEL.output_bed,
            svtype = "DEL",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call Helpers.BedtoolsClosest as CompareDUP {
        input:
            bed_a = SplitEval.dup_bed,
            bed_b = SplitTruth.dup_bed,
            svtype = "DUP",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }
    
    call Helpers.SelectMatchedSVs as CalcuDUP {
        input:
            input_bed = CompareDUP.output_bed,
            svtype = "DUP",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call Helpers.BedtoolsClosest as CompareINS {
        input:
            bed_a = SplitEval.ins_bed,
            bed_b = SplitTruth.ins_bed,
            svtype = "INS",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }
    
    call Helpers.SelectMatchedINSs as CalcuINS {
        input:
            input_bed = CompareINS.output_bed,
            svtype = "INS",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call Helpers.BedtoolsClosest as CompareINV {
        input:
            bed_a = SplitEval.inv_bed,
            bed_b = SplitTruth.inv_bed,
            svtype = "INV",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }
    
    call Helpers.SelectMatchedSVs as CalcuINV {
        input:
            input_bed = CompareINV.output_bed,
            svtype = "INV",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call Helpers.BedtoolsClosest as CompareBND {
        input:
            bed_a = SplitEval.bnd_bed,
            bed_b = SplitTruth.bnd_bed,
            svtype = "BND",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }
    
    call Helpers.SelectMatchedINSs as CalcuBND {
        input:
            input_bed = CompareBND.output_bed,
            svtype = "BND",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call Helpers.ConcatFiles as MergeComparisons {
        input:
            files = [CalcuDEL.output_comp, CalcuDUP.output_comp, CalcuINS.output_comp, CalcuINV.output_comp, CalcuBND.output_comp],
            outfile_name = "~{prefix}.comparison.bed",
            docker_image = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File closest_bed = MergeComparisons.concatenated_file
    }
} 