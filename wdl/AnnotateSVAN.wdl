version 1.0
    
import "Structs.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_index
        String prefix
        
        File vntr_bed
        File exons_bed
        File repeats_bed
        File consensus_fasta
        File reference_fasta
        File reference_fasta_fai
        
        String svan_docker
        
        RuntimeAttr? runtime_attr_generate_trf
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_merge
    }

    call SeparateInsertionsDeletions {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call GenerateTRFForInsertions {
        input:
            vcf = SeparateInsertionsDeletions.ins_vcf,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_generate_trf
    }

    call AnnotateInsertions {
        input:
            vcf = SeparateInsertionsDeletions.ins_vcf,
            trf_output = GenerateTRFForInsertions.trf_output,
            vntr_bed = vntr_bed,
            exons_bed = exons_bed,
            repeats_bed = repeats_bed,
            consensus_fasta = consensus_fasta,
            reference_fasta = reference_fasta,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_annotate_ins
    }

    call GenerateTRFForDeletions {
        input:
            vcf = SeparateInsertionsDeletions.del_vcf,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_generate_trf
    }

    call AnnotateDeletions {
        input:
            vcf = SeparateInsertionsDeletions.del_vcf,
            trf_output = GenerateTRFForDeletions.trf_output,
            vntr_bed = vntr_bed,
            exons_bed = exons_bed,
            repeats_bed = repeats_bed,
            consensus_fasta = consensus_fasta,
            reference_fasta = reference_fasta,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_annotate_del
    }

    call MergeAnnotatedVcfs {
        input:
            original_vcf = vcf,
            annotated_ins_vcf = AnnotateInsertions.annotated_vcf,
            annotated_del_vcf = AnnotateDeletions.annotated_vcf,
            other_variants_vcf = SeparateInsertionsDeletions.other_vcf,
            prefix = prefix,
            svan_docker = svan_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File annotated_vcf = MergeAnnotatedVcfs.merged_vcf
        File annotated_vcf_index = MergeAnnotatedVcfs.merged_vcf_index
        
        File? ins_vcf = SeparateInsertionsDeletions.ins_vcf
        File? del_vcf = SeparateInsertionsDeletions.del_vcf
        File other_vcf = SeparateInsertionsDeletions.other_vcf
        
        File? ins_trf_output = GenerateTRFForInsertions.trf_output
        File? del_trf_output = GenerateTRFForDeletions.trf_output
        
        File? annotated_ins_vcf = AnnotateInsertions.annotated_vcf
        File? annotated_del_vcf = AnnotateDeletions.annotated_vcf
    }
}

task SeparateInsertionsDeletions {
    input {
        File vcf
        File vcf_index
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include 'INFO/SVTYPE="INS" && INFO/SVLEN > 0' \
            -O z -o ~{prefix}.insertions.vcf.gz
        
        bcftools view ~{vcf} \
            --include 'INFO/SVTYPE="DEL" && INFO/SVLEN < 0' \
            -O z -o ~{prefix}.deletions.vcf.gz
        
        bcftools view ~{vcf} \
            --exclude '(INFO/SVTYPE="INS" && INFO/SVLEN > 0) || (INFO/SVTYPE="DEL" && INFO/SVLEN < 0)' \
            -O z -o ~{prefix}.other_variants.vcf.gz
        
        tabix -p vcf ~{prefix}.insertions.vcf.gz
        tabix -p vcf ~{prefix}.deletions.vcf.gz
        tabix -p vcf ~{prefix}.other_variants.vcf.gz
    >>>

    output {
        File ins_vcf = "~{prefix}.insertions.vcf.gz"
        File ins_vcf_index = "~{prefix}.insertions.vcf.gz.tbi"
        File del_vcf = "~{prefix}.deletions.vcf.gz"
        File del_vcf_index = "~{prefix}.deletions.vcf.gz.tbi"
        File other_vcf = "~{prefix}.other_variants.vcf.gz"
        File other_vcf_index = "~{prefix}.other_variants.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([vcf, vcf_index], "GB") * 4) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GenerateTRFForInsertions {
    input {
        File vcf
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        python3 /app/SVAN/scripts/ins2fasta.py ~{vcf} work_dir
        
        trf work_dir/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.ins_trf.out
    >>>

    output {
        File trf_output = "~{prefix}.ins_trf.out"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB") * 10) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GenerateTRFForDeletions {
    input {
        File vcf
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        python3 /app/SVAN/scripts/del2fasta.py ~{vcf} work_dir
        
        trf work_dir/deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.del_trf.out
    >>>

    output {
        File trf_output = "~{prefix}.del_trf.out"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB") * 10) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateInsertions {
    input {
        File vcf
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File consensus_fasta
        File reference_fasta
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        python3 /app/SVAN/SVAN-INS.py \
            ~{vcf} \
            ~{trf_output} \
            ~{vntr_bed} \
            ~{exons_bed} \
            ~{repeats_bed} \
            ~{consensus_fasta} \
            ~{reference_fasta} \
            ~{prefix}.svan_annotated \
            -o work_dir

        bgzip -c work_dir/~{prefix}.svan_annotated.vcf > ~{prefix}.ins_annotated.vcf.gz
        tabix -p vcf ~{prefix}.ins_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.ins_annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.ins_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: ceil(size([vcf, reference_fasta, consensus_fasta], "GB") * 8) + 100,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateDeletions {
    input {
        File vcf
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File consensus_fasta
        File reference_fasta
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        python3 /app/SVAN/SVAN-DEL.py \
            ~{vcf} \
            ~{trf_output} \
            ~{vntr_bed} \
            ~{exons_bed} \
            ~{repeats_bed} \
            ~{consensus_fasta} \
            ~{reference_fasta} \
            ~{prefix}.svan_annotated \
            -o work_dir

        bgzip -c work_dir/~{prefix}.svan_annotated.vcf > ~{prefix}.del_annotated.vcf.gz
        tabix -p vcf ~{prefix}.del_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.del_annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.del_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: ceil(size([vcf, reference_fasta, consensus_fasta], "GB") * 8) + 100,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeAnnotatedVcfs {
    input {
        File original_vcf
        File? annotated_ins_vcf
        File? annotated_del_vcf
        File other_variants_vcf
        String prefix
        String svan_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Create list of VCFs to merge
        echo "~{other_variants_vcf}" > vcf_list.txt
        
        if [[ "~{annotated_ins_vcf}" != "" ]]; then
            echo "~{annotated_ins_vcf}" >> vcf_list.txt
        fi
        
        if [[ "~{annotated_del_vcf}" != "" ]]; then
            echo "~{annotated_del_vcf}" >> vcf_list.txt
        fi

        # Sort list for consistent merging
        sort vcf_list.txt > sorted_vcf_list.txt
        
        # Merge VCFs maintaining original header and variant order
        bcftools concat -f sorted_vcf_list.txt -O z -o ~{prefix}.svan_annotated.vcf.gz
        
        tabix -p vcf ~{prefix}.svan_annotated.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.svan_annotated.vcf.gz"
        File merged_vcf_index = "~{prefix}.svan_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size([original_vcf, annotated_ins_vcf, annotated_del_vcf, other_variants_vcf], "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: svan_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}