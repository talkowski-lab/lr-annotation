version 1.0

import "general/Helpers.wdl"

workflow TRGT {
    input {
        File bam
        File bai

        String prefix
        String sample_id
        String sex

        File ref_fa
        File ref_fai
        File repeatCatalog
        String catalog_name

        String gcs_out_dir
        String trgt_docker
        RuntimeAttr? runtime_attr_process_with_trgt
    }

    Boolean is_female = 'F' == sex
    String outdir = sub(gcs_out_dir, "/$", "") + "/TRGT/" + catalog_name

    call ProcessWithTRGT {
        input:
            bam = bam,
            bai = bai,
            is_female = is_female,
            outprefix = sample_id,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            repeatCatalog = repeatCatalog,
            catalog_name = catalog_name,
            docker = trgt_docker,
            runtime_attr_override = runtime_attr_process_with_trgt
    }

    call Helpers.FinalizeToFile as FinalizeTrgtVCF { 
        input: 
            outdir = outdir, 
            file = ProcessWithTRGT.trgt_output_vcf     
    }

    call Helpers.FinalizeToFile as FinalizeTrgtTBI { 
        input: 
            outdir = outdir, 
            file = ProcessWithTRGT.trgt_output_vcf_idx 
    }

    output {
        File trgt_vcf = FinalizeTrgtVCF.gcs_path
        File trgt_vcf_idx = FinalizeTrgtTBI.gcs_path
        # File trgt_bam = ProcessWithTRGT.trgt_output_bam
    }
}

task ProcessWithTRGT {
    input {
        File bam
        File bai
        String outprefix
        File repeatCatalog
        String catalog_name
        Boolean is_female
        Boolean verbose = false
        File ref_fa
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String vcf_out_name = outprefix + "_trgt." + catalog_name
    String karyotype = if(is_female) then "XX" else "XY"

    command <<<
        set -euo pipefail

        nproc=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        time \
        trgt \
            ~{true='--verbose' false=' ' verbose} \
            genotype \
            --genome ~{ref_fa} \
            --repeats ~{repeatCatalog} \
            --reads ~{bam} \
            --threads "${nproc}" \
            --output-prefix ~{vcf_out_name} \
            --karyotype ~{karyotype}

        find . | sed -e "s/[^-][^\/]*\// |/g" -e "s/|\([^ ]\)/|-\1/"

        bcftools sort \
            -Ob \
            -o ~{vcf_out_name}.sorted.vcf.gz \
            ~{vcf_out_name}.vcf.gz
        
        mv ~{vcf_out_name}.sorted.vcf.gz \
            ~{vcf_out_name}.vcf.gz
        
        bcftools index \
            -t \
            ~{vcf_out_name}.vcf.gz
    >>>

    output {
        File trgt_output_vcf = "~{vcf_out_name}.vcf.gz"
        File trgt_output_vcf_idx = "~{vcf_out_name}.vcf.gz.tbi"
        # File trgt_output_bam = "~{vcf_out_name}.spanning.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 6,
        disk_gb: 10 + 1.3 * ceil(size(bam, "GiB") + size(ref_fa, "GiB") + size(repeatCatalog, "GiB")),
        boot_disk_gb: 20,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
