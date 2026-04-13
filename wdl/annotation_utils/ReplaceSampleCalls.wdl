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

sample_variant_ids = set()
with pysam.VariantFile("~{sample_vcf}") as s:
    sample_name = list(s.header.samples)[0]
    for record in s:
        if record.id:
            sample_variant_ids.add(record.id)

vcf_in = pysam.VariantFile("~{cohort_vcf}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=vcf_in.header)

for record in vcf_in:
    allele_type = record.info.get('allele_type', None)
    if allele_type == 'trv' and record.id:
        sample = record.samples[sample_name]
        gt = sample['GT']
        is_hom_ref = gt is not None and len(gt) > 0 and all(a == 0 for a in gt)
        if not is_hom_ref:
            id_prefix = record.id.split('_')[0]
            if id_prefix not in sample_variant_ids:
                sample['GT'] = tuple(None for _ in gt)
                if 'PS' in sample:
                    sample['PS'] = None
    vcf_out.write(record)

vcf_in.close()
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
