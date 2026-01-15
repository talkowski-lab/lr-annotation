version 1.0

import "../utils/Structs.wdl"

workflow SubsetVcfToSamples {
    input {
        File vcf
        File vcf_idx
        Array[String] samples

        String prefix
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
    }

    call SubsetVcfToSampleList {
        input:
            vcf = vcf,
            vcf_index = vcf_idx,
            samples = samples,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    output {
        File subset_vcf = SubsetVcfToSampleList.subset_vcf
        File subset_vcf_idx = SubsetVcfToSampleList.subset_vcf_index
    }
}

task SubsetVcfToSampleList {
    input {
        File vcf
        File vcf_index
        Array[String] samples
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > samples.txt <<EOF
~{sep='\n' samples}
EOF

        bcftools view \
            --samples-file samples.txt \
            --min-ac 1 \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
