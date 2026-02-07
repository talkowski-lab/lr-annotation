version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TRGTLPS {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        
        String trgt_lps_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_trgt_lps
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_finalize
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call RunTRGTLPS {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.{contig}.lps",
                docker = trgt_lps_docker,
                runtime_attr_override = runtime_attr_trgt_lps
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = RunTRGTLPS.lps_tsv,
            prefix = "~{prefix}.concat",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call AddHeader {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            tsv = ConcatTsvs.concatenated_tsv,
            prefix =  "~{prefix}.headered",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_finalize
    }

    output {
        File trgt_lps_tsv = AddHeader.tsv_with_header
    }
}

task RunTRGTLPS {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    command <<<
        set -eou pipefail
        
        trgt-lps \
            --vcf ~{vcf} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            > ~{prefix}.tsv
    >>>

    output {
        File lps_tsv = "~{prefix}.tsv"
    }
}

task AddHeader {
    input {
        File vcf
        File vcf_idx
        File tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail
        
        bcftools query -l ~{vcf} > samples.txt
        
        echo -e "ID\tmotif\t$(cat samples.txt | tr '\n' '\t' | sed 's/\t$//')" > header.txt
        
        cat header.txt ~{tsv} > ~{prefix}.lps.tsv
    >>>

    output {
        File tsv_with_header = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(vcf, "GB")) + 2 * ceil(size(tsv, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
