version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_truth
        String prefix
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_concat
    }

    call Helpers.SplitQueryVcf as SplitEval {
        input:
            vcf = vcf_eval,
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

    # ... similar calls for DUP, INS, INV, BND ...

    call Helpers.ConcatFiles as MergeBeds {
        input:
            files = [CompareDEL.output_bed, ...],
            outfile_name = "~{prefix}.closest.bed",
            docker_image = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File closest_bed = MergeBeds.concatenated_file
    }
} 