version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateDuplications {
    input {
        File vcf
        File vcf_idx
        File annotations_tsv
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetTsvToContig {
            input:
                tsv = annotations_tsv,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tsv
        }

        call UpdateAlleleType {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                annotations_tsv = SubsetTsvToContig.subset_tsv,
                prefix = "~{prefix}.~{contig}.dup_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = UpdateAlleleType.annotated_vcf,
            vcfs_idx = UpdateAlleleType.annotated_vcf_idx,
            prefix = prefix + ".dup_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File dup_annotated_vcf = ConcatVcfs.concat_vcf
        File dup_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task UpdateAlleleType {
    input {
        File vcf
        File vcf_idx
        File annotations_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bgzip -c "~{annotations_tsv}" > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        
        bcftools annotate \
            -a annotations.tsv.gz \
            -c CHROM,POS,REF,ALT,ID,INFO/SVTYPE \
            -Oz -o ~{prefix}.annotated.vcf.gz \
            "~{vcf}"

        tabix -p vcf ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(annotations_tsv, "GB")) + 20,
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
