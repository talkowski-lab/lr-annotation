version 1.0

import "general/Structs.wdl"

workflow AnnotateRepeatMasker {
    input {
        File vcf
        File vcf_idx

        String prefix

        String utils_docker
        String repeatmasker_docker
        RuntimeAttr? runtime_attr_vcf_ins_to_fa
        RuntimeAttr? runtime_attr_repeat_masker
    }

    call VCF_INS_to_fa {
        input:
            vcf = vcf,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_vcf_ins_to_fa
    }

    call RepeatMasker {
        input:
            fasta = VCF_INS_to_fa.INS_fa,
            prefix = prefix,
            docker = repeatmasker_docker,
            runtime_attr_override = runtime_attr_repeat_masker
    }

    output {
        File RM_out = RepeatMasker.RMout
        File RM_fa = VCF_INS_to_fa.INS_fa
    }
}

task VCF_INS_to_fa {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -i 'INFO/SVLEN>=50 && INFO/SVTYPE=="INS"' ~{vcf} | bcftools view -e 'ALT ~ "<"' > ~{prefix}.INS.vcf

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ~{prefix}.INS.vcf | awk 'length($3)==1 {print ">"$1":"$2";"$3"\n"$4}'> ~{prefix}_INS.tmp.fa
        
        seqkit rename -N1 ~{prefix}_INS.tmp.fa > ~{prefix}_INS.fa
    >>>

    output {
        File INS_fa = '~{prefix}_INS.fa'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2*ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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

task RepeatMasker {
    input {
        File fasta
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
  
        RepeatMasker \
            -e rmblast \
            -species human \
            -pa 4 \
            -s \
            ~{fasta}
        
        mv ~{fasta}.out ./~{prefix}.out
    >>>

    output {
        File RMout = '~{prefix}.out'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 5*ceil(size(fasta, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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
