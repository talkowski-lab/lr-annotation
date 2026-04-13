version 1.0

import "../utils/Structs.wdl"

workflow ReplaceSampleCalls {
    input {
        File sample_vcf
        File sample_vcf_idx
        File cohort_vcf
        File cohort_vcf_idx
        String prefix

        String docker

        RuntimeAttr? runtime_attr_replace
    }

    call ReplaceCalls {
        input:
            sample_vcf = sample_vcf,
            sample_vcf_idx = sample_vcf_idx,
            cohort_vcf = cohort_vcf,
            cohort_vcf_idx = cohort_vcf_idx,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_replace
    }

    output {
        File replaced_vcf = ReplaceCalls.replaced_vcf
        File replaced_vcf_idx = ReplaceCalls.replaced_vcf_idx
    }
}

task ReplaceCalls {
    input {
        File sample_vcf
        File sample_vcf_idx
        File cohort_vcf
        File cohort_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

sample_in = pysam.VariantFile("~{sample_vcf}")
sample_name = list(sample_in.header.samples)[0]

sample_variant_ids = {}
for record in sample_in:
    data = { 'GT': record.samples[sample_name]['GT'] }
    if 'PS' in record.format:
        data['PS'] = record.samples[sample_name]['PS']
    if 'PF' in record.format:
        data['PF'] = record.samples[sample_name]['PF']
    sample_variant_ids[record.id] = data
sample_in.close()

cohort_in = pysam.VariantFile("~{cohort_vcf}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=cohort_in.header)

for record in cohort_in:
    match_data = sample_variant_ids.get(record.id.split('_')[0])

    if match_data is not None:
        # Set GT/PS/PF from sample_vcf for matched variants
        record.samples[sample_name]['GT'] = match_data['GT']
        if 'PS' in match_data:
            record.samples[sample_name]['PS'] = match_data['PS']
        if 'PF' in match_data:
            record.samples[sample_name]['PF'] = match_data['PF']
    elif record.info.get('allele_type') == 'trv':
        # Clear GT/PS/PF for unmatched trv variants
        record.samples[sample_name]['GT'] = tuple(None for _ in record.samples[sample_name]['GT'])
        if 'PS' in record.format:
            record.samples[sample_name]['PS'] = None
        if 'PF' in record.format:
            record.samples[sample_name]['PF'] = None

    vcf_out.write(record)

cohort_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File replaced_vcf = "~{prefix}.vcf.gz"
        File replaced_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(cohort_vcf, "GB") + size(sample_vcf, "GB")) + 20,
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
