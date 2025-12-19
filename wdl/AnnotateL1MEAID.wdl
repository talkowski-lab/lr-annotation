version 1.0

import "general/Structs.wdl"
import "RepeatMasker.wdl" as RM

workflow AnnotateL1MEAID {
    input {
        File vcf
        File vcf_idx
        String prefix

        String utils_docker
        String repeatmasker_docker
        String annotate_l1meaid_docker
        String annotate_l1meaid_filter_docker

        RuntimeAttr? runtime_attr_ins_to_fa
        RuntimeAttr? runtime_attr_repeat_masker
        RuntimeAttr? runtime_attr_limeaid
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_annotate
    }

    call RM.RepeatMasker {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            utils_docker = utils_docker,
            repeatmasker_docker = repeatmasker_docker,
            runtime_attr_ins_to_fa = runtime_attr_ins_to_fa,
            runtime_attr_repeat_masker = runtime_attr_repeat_masker
    }

    call L1MEAID {
        input:
            rm_fa = RepeatMasker.rm_fa,
            rm_out = RepeatMasker.rm_out,
            prefix = prefix,
            docker = annotate_l1meaid_docker,
            runtime_attr_override = runtime_attr_limeaid
    }

    call L1MEAIDFilter {
        input:
            limeaid_output = L1MEAID.limeaid_output,
            prefix = prefix,
            docker = annotate_l1meaid_filter_docker,
            runtime_attr_override = runtime_attr_filter
    }

    call AnnotateVCF {
        input:
            vcf = vcf,
            l1meaid_filtered_tsv = L1MEAIDFilter.filtered_output,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_annotate
    }

    output {
        File annotated_vcf = AnnotateVCF.annotated_vcf
        File annotated_vcf_idx = AnnotateVCF.annotated_vcf_idx
        File l1meaid_filtered_output = L1MEAIDFilter.filtered_output
        File rm_out = RepeatMasker.rm_out
    }
}

task L1MEAID {
    input {
        File rm_fa
        File rm_out
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/src/L1ME-AID/limeaid.py \
            -i ~{rm_fa} \
            -r ~{rm_out} \
            -o ~{prefix}_limeaid_output.txt
    >>>

    output {
        File limeaid_output = "~{prefix}_limeaid_output.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2*ceil(size(rm_fa, "GB")) + 5,
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

task L1MEAIDFilter {
    input {
        File limeaid_output
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        perl /opt/src/utility/limeaid.filter.pl \
            ~{limeaid_output} \
            > ~{prefix}_filtered.tsv
    >>>

    output {
        File filtered_output = "~{prefix}_filtered.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2*ceil(size(limeaid_output, "GB")) + 5,
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

task AnnotateVCF {
    input {
        File vcf
        File l1meaid_filtered_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
import sys

input_tsv = "~{l1meaid_filtered_tsv}"
output_anno = "annotations.tsv"

with open(input_tsv, 'r') as f_in, open(output_anno, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) < 12: continue
        
        var_id_str = parts[0]
        sequence = parts[1]
        
        classification = parts[8]
        structure = parts[10]
        
        me_type = None
        
        if classification == "SINE/Alu" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "ALU"
        elif classification == "Retroposon/SVA" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "SVA"
        elif classification == "LINE/L1" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "LINE"
        
        if me_type:
            id_parts = var_id_str.split(';')
            if len(id_parts) >= 2:
                location_part = id_parts[0]
                ref = id_parts[1]
                
                try:
                    chrom, pos = location_part.rsplit(':', 1)
                    f_out.write(f"{chrom}\t{pos}\t{ref}\t{sequence}\t{me_type}\n")
                except ValueError:
                    sys.stderr.write(f"Skipping malformed location: {var_id_str}\n")
EOF
        
        if bcftools view -h ~{vcf} | grep -q "ID=ME_TYPE"; then
            cp ~{vcf} headered.vcf.gz
        else
            echo '##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">' > header.txt
            bcftools annotate -h header.txt -Oz -o headered.vcf.gz ~{vcf}
        fi
        bcftools index -t headered.vcf.gz
        
        if [ -s annotations.tsv ]; then
            sort -k1,1 -k2,2n annotations.tsv | bgzip -c > annotations.tsv.gz
            tabix -s1 -b2 -e2 annotations.tsv.gz
            
            bcftools annotate \
                -a annotations.tsv.gz \
                -c CHROM,POS,REF,ALT,INFO/ME_TYPE \
                -Oz -o ~{prefix}.l1meaid_annotated.vcf.gz \
                headered.vcf.gz
            
            bcftools index -t ~{prefix}.l1meaid_annotated.vcf.gz
        else
            echo "No matching L1MEAID variants found."
            mv headered.vcf.gz ~{prefix}.l1meaid_annotated.vcf.gz
            bcftools index -t ~{prefix}.l1meaid_annotated.vcf.gz
        fi
    >>>

    output {
        File annotated_vcf = "~{prefix}.l1meaid_annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.l1meaid_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 5,
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
