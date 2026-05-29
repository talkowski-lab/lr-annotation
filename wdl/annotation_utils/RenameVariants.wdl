version 1.0

import "../utils/Structs.wdl"

workflow RenameVariants {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker

        RuntimeAttr? runtime_attr_rename
    }

    call RenameVariantIds {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_rename
    }

    output {
        File renamed_vcf = RenameVariantIds.renamed_vcf
        File renamed_vcf_idx = RenameVariantIds.renamed_vcf_idx
    }
}

task RenameVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import re
import pysam

def clean(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value))

vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("~{prefix}.unsorted.vcf.gz", "w", header=vcf_in.header)

counter = 0
for record in vcf_in:
    counter += 1
    svtype = record.info.get("SVTYPE", "NA")
    if isinstance(svtype, (list, tuple)):
        svtype = svtype[0] if svtype else "NA"
    record.id = f"{clean(record.chrom)}_{record.pos}_{clean(svtype)}_{counter}"
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        bcftools sort \
            -Oz -o ~{prefix}.vcf.gz \
            ~{prefix}.unsorted.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
