version 1.0

import "../utils/Structs.wdl"

workflow SplitVcfPerContig {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Boolean create_no_geno = false
        Boolean modify_snv_ids = false

        String utils_docker

        RuntimeAttr? runtime_attr_split_vcf
    }

    call SplitByContig {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            contigs = contigs,
            prefix = prefix,
            create_no_geno = create_no_geno,
            modify_snv_ids = modify_snv_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    output {
        Array[File] contig_vcfs = SplitByContig.contig_vcfs
        Array[File] contig_vcf_idxs = SplitByContig.contig_vcf_idxs
        Array[File] contig_no_geno_vcfs = SplitByContig.contig_no_geno_vcfs
        Array[File] contig_no_geno_vcf_idxs = SplitByContig.contig_no_geno_vcf_idxs
    }
}

task SplitByContig {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        Boolean create_no_geno
        Boolean modify_snv_ids
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch contig_vcfs.list contig_vcf_idxs.list contig_no_geno_vcfs.list contig_no_geno_vcf_idxs.list

        while IFS= read -r contig; do
            if ~{modify_snv_ids}; then
                bcftools view \
                    -r "${contig}" \
                    ~{vcf} \
                | awk 'BEGIN{OFS="\t"} /^#/ {print; next} $8 ~ /(^|;)allele_type=snv(;|$)/ {$3=$1"-"$2"-"$4"-"$5} {print}' \
                | bgzip -c > "~{prefix}.${contig}.vcf.gz"
            else
                bcftools view \
                    -r "${contig}" \
                    ~{vcf} \
                    -Oz -o "~{prefix}.${contig}.vcf.gz"
            fi

            tabix -p vcf "~{prefix}.${contig}.vcf.gz"

            echo "~{prefix}.${contig}.vcf.gz" >> contig_vcfs.list
            echo "~{prefix}.${contig}.vcf.gz.tbi" >> contig_vcf_idxs.list

            if ~{create_no_geno}; then
                bcftools view \
                    -G \
                    "~{prefix}.${contig}.vcf.gz" \
                    -Oz -o "~{prefix}.${contig}.no_geno.vcf.gz"
                
                tabix -p vcf "~{prefix}.${contig}.no_geno.vcf.gz"

                echo "~{prefix}.${contig}.no_geno.vcf.gz" >> contig_no_geno_vcfs.list
                echo "~{prefix}.${contig}.no_geno.vcf.gz.tbi" >> contig_no_geno_vcf_idxs.list
            fi
        done < ~{write_lines(contigs)}
    >>>

    output {
        Array[File] contig_vcfs = read_lines("contig_vcfs.list")
        Array[File] contig_vcf_idxs = read_lines("contig_vcf_idxs.list")
        Array[File] contig_no_geno_vcfs = read_lines("contig_no_geno_vcfs.list")
        Array[File] contig_no_geno_vcf_idxs = read_lines("contig_no_geno_vcf_idxs.list")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
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
