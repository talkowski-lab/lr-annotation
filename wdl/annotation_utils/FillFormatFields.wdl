version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        String prefix

        Array[String] format_fields
        Boolean fill_sv_gts = false
        String? sv_source_tag

        String utils_docker

        RuntimeAttr? runtime_attr_fill
    }

    call Helpers.FillVcfFormatFields {
        input:
            unfilled_vcf = unfilled_vcf,
            unfilled_vcf_idx = unfilled_vcf_idx,
            filled_vcf = filled_vcf,
            filled_vcf_idx = filled_vcf_idx,
            format_fields = format_fields,
            fill_sv_gts = fill_sv_gts,
            sv_source_tag = sv_source_tag,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_fill
    }

    output {
        File output_vcf = FillVcfFormatFields.output_vcf
        File output_vcf_idx = FillVcfFormatFields.output_vcf_idx
    }
}
