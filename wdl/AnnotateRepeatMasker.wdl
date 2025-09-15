version 1.0

import "Structs.wdl"

workflow AnnotateRepeatMasker {
    input {
        File vcf
        File vcf_index

        RuntimeAttr? runtime_attr_override_vcf_ins_to_fa
        RuntimeAttr? runtime_attr_override_repeat_masker
    }

    call VCF_INS_to_fa { 
        input: 
            vcf = vcf,
            runtime_attr_override = runtime_attr_override_vcf_ins_to_fa
    }

    call RepeatMasker {
        input: 
            fasta = VCF_INS_to_fa.INS_fa,
            runtime_attr_override = runtime_attr_override_repeat_masker
    }

    output {
        File RMout = RepeatMasker.RMout
        File INSfa = VCF_INS_to_fa.INS_fa
    }
}

task VCF_INS_to_fa {
    input {
        File vcf
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcf, "GB"))+20
    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euxo pipefail

        bcftools view -i 'INFO/SVLEN>=50 && INFO/SVTYPE=="INS"' ~{vcf} | bcftools view -e 'ALT ~ "<"' > ~{prefix}.INS.vcf
        # the next line filters out "INS" calls in the integrated callset that are actually INV
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ~{prefix}.INS.vcf |awk 'length($3)==1 {print ">"$1":"$2";"$3"\n"$4}'> ~{prefix}_INS.tmp.fa
        
        # rename duplicate fasta IDs
        seqkit rename -N1 ~{prefix}_INS.tmp.fa > ~{prefix}_INS.fa
    >>>

    output {
        File INS_fa = '~{prefix}_INS.fa'
    }

     RuntimeAttr default_attr = object {
         cpu_cores:          1,
         mem_gb:             2,
         disk_gb:            disk_size,
         boot_disk_gb:       10,
         preemptible_tries:  3,
         max_retries:        0,
         docker:             "quay.io/ymostovoy/lr-utils-basic:latest"
     }

     RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
     runtime {
         cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
         memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
         disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
         bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
         preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
         maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
         docker:                 select_first([runtime_attr.docker,            default_attr.docker])
     }
}

task RepeatMasker {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 5*ceil(size(fasta, "GB"))+20
    String prefix = basename(fasta)

    command <<<
        set -euxo pipefail
  
        RepeatMasker -e rmblast -pa 4 -s -species human ~{fasta}
        mv ~{fasta}.out ./~{prefix}.out
     >>>

     output {
        File RMout = '~{prefix}.out'
     }

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "dfam/tetools:1.8"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
