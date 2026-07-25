version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterTRGTCalls {
    input {
        File trgt_vcf
        File trgt_vcf_idx
        String prefix

        Int? min_repeat_unit
        Int? min_length_diff
        Int? max_catalog_length

        String utils_docker

        RuntimeAttr? runtime_attr_filter
    }

    call Helpers.FilterTRGTVcf {
        input:
            vcf = trgt_vcf,
            vcf_idx = trgt_vcf_idx,
            min_repeat_unit = min_repeat_unit,
            min_length_diff = min_length_diff,
            max_catalog_length = max_catalog_length,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File trgt_filtered_vcf = FilterTRGTVcf.processed_vcf
        File trgt_filtered_vcf_idx = FilterTRGTVcf.processed_vcf_idx
    }
}
