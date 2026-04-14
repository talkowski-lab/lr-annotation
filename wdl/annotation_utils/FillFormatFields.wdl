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
        String? include_field
        String? include_value
        Boolean fill_alt_gts = false
        Boolean fill_ref_gts = false

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
            include_field = include_field,
            include_value = include_value,
            fill_alt_gts = fill_alt_gts,
            fill_ref_gts = fill_ref_gts,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_fill
    }

    output {
        File refilled_vcf = FillVcfFormatFields.output_vcf
        File refilled_vcf_idx = FillVcfFormatFields.output_vcf_idx
    }
}
