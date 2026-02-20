version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateVRS {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        
        File seqrepo_tar

        String utils_docker
        String vrs_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_annotate_vrs
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call AnnotateVcfWithVRS {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                seqrepo_tar_gz = seqrepo_tar_gz,
                prefix = "~{prefix}.~{contig}.vrs",
                docker = vrs_docker,
                runtime_attr_override = runtime_attr_annotate_vrs
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateVcfWithVRS.annotated_vcf,
            vcf_idxs = AnnotateVcfWithVRS.annotated_vcf_idx,
            prefix = "~{prefix}.vrs_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File vrs_annotated_vcf = ConcatVcfs.concat_vcf
        File vrs_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateVcfWithVRS {
    input {
        File vcf
        File vcf_idx
        File seqrepo_tar_gz
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Setup SeqRepo
        mkdir seqrepo_data
        tar -xzf ~{seqrepo_tar_gz} -C seqrepo_data --strip-components 1
        
        # Locate the directory containing the aliases.sqlite3 file to build the URI
        SEQREPO_PATH=$(find $(pwd)/seqrepo_data -name "aliases.sqlite3" | xargs dirname)
        echo "SeqRepo path determined: $SEQREPO_PATH"

        # Run vrs-annotate
        vrs-annotate vcf \
            --dataproxy-uri="seqrepo+file://${SEQREPO_PATH}" \
            --vcf-out ~{prefix}.vcf.gz \
            --vrs-attributes \
            ~{vcf}

        # Index output
        tabix -p vcf ~{prefix}.vcf.gz
        
        # Clean up
        rm -rf seqrepo_data
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4, 
        disk_gb: 2 * ceil(size(vcf, "GB") + size(seqrepo_tar_gz, "GB")) + 20,
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
