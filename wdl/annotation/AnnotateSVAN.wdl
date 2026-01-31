version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        
        File vntr_bed
        File exons_bed
        File repeats_bed

        File mei_fa
        File mei_fa_amb
        File mei_fa_ann
        File mei_fa_bwt
        File mei_fa_pac
        File mei_fa_sa
        File mei_fa_mmi
        File ref_fa

        String annotate_svan_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_separate
        RuntimeAttr? runtime_attr_generate_trf_ins
        RuntimeAttr? runtime_attr_generate_trf_del
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_concat_ins
        RuntimeAttr? runtime_attr_concat_del
        RuntimeAttr? runtime_attr_concat_final
        RuntimeAttr? runtime_attr_merge_headers
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                drop_genotypes = true,
                prefix = prefix,
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ResetVcfFilters {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.reset_filters",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call SeparateInsertionsDeletions {
            input:
                vcf = ResetVcfFilters.reset_vcf,
                vcf_idx = ResetVcfFilters.reset_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_separate
        }

        # Insertions
        call GenerateTRF as GenerateTRFForInsertions {
            input:
                vcf = SeparateInsertionsDeletions.ins_vcf,
                vcf_idx = SeparateInsertionsDeletions.ins_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                mode = "ins",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_generate_trf_ins
        }

        call RunSvanAnnotation as AnnotateInsertions {
            input:
                vcf = SeparateInsertionsDeletions.ins_vcf,
                vcf_idx = SeparateInsertionsDeletions.ins_vcf_idx,
                trf_output = GenerateTRFForInsertions.trf_output,
                vntr_bed = vntr_bed,
                exons_bed = exons_bed,
                repeats_bed = repeats_bed,
                mei_fa = mei_fa,
                mei_fa_amb = mei_fa_amb,
                mei_fa_ann = mei_fa_ann,
                mei_fa_bwt = mei_fa_bwt,
                mei_fa_pac = mei_fa_pac,
                mei_fa_sa = mei_fa_sa,
                mei_fa_mmi = mei_fa_mmi,
                ref_fa = ref_fa,
                prefix = "~{prefix}.~{contig}",
                mode = "ins",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_annotate_ins
        }

        call Helpers.ExtractVcfAnnotations as ExtractIns {
            input:
                vcf = AnnotateInsertions.annotated_vcf,
                vcf_idx = AnnotateInsertions.annotated_vcf_idx,
                original_vcf = SeparateInsertionsDeletions.ins_vcf,
                original_vcf_idx = SeparateInsertionsDeletions.ins_vcf_idx,
                add_header_row = true,
                prefix = "~{prefix}.~{contig}.ins",
                docker = annotate_svan_docker
        }

        # Deletions
        call GenerateTRF as GenerateTRFForDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_idx = SeparateInsertionsDeletions.del_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                mode = "del",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_generate_trf_del
        }

        call RunSvanAnnotation as AnnotateDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_idx = SeparateInsertionsDeletions.del_vcf_idx,
                trf_output = GenerateTRFForDeletions.trf_output,
                vntr_bed = vntr_bed,
                exons_bed = exons_bed,
                repeats_bed = repeats_bed,
                mei_fa = mei_fa,
                mei_fa_amb = mei_fa_amb,
                mei_fa_ann = mei_fa_ann,
                mei_fa_bwt = mei_fa_bwt,
                mei_fa_pac = mei_fa_pac,
                mei_fa_sa = mei_fa_sa,
                mei_fa_mmi = mei_fa_mmi,
                ref_fa = ref_fa,
                prefix = "~{prefix}.~{contig}",
                mode = "del",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_annotate_del
        }

        call Helpers.ExtractVcfAnnotations as ExtractDel {
            input:
                vcf = AnnotateDeletions.annotated_vcf,
                vcf_idx = AnnotateDeletions.annotated_vcf_idx,
                original_vcf = SeparateInsertionsDeletions.del_vcf,
                original_vcf_idx = SeparateInsertionsDeletions.del_vcf_idx,
                add_header_row = true,
                prefix = "~{prefix}.~{contig}.del",
                docker = annotate_svan_docker
        }
    }
    
    call Helpers.ConcatAlignedTsvs {
        input:
            tsvs = flatten([ExtractIns.annotations_tsv, ExtractDel.annotations_tsv]),
            prefix = prefix + ".svan_annotations",
            docker = annotate_svan_docker,
            runtime_attr_override = runtime_attr_concat_final
    }

    output {
        File annotations_tsv_svan = ConcatAlignedTsvs.merged_tsv
        File annotations_header_svan = ConcatAlignedTsvs.merged_header
    }
}

task SeparateInsertionsDeletions {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include 'INFO/allele_type="INS"' \
            -Oz -o ~{prefix}.insertions.vcf.gz
        
        bcftools view ~{vcf} \
            --include 'INFO/allele_type="DEL"' \
            -Oz -o ~{prefix}.deletions.vcf.gz
        
        tabix -p vcf ~{prefix}.insertions.vcf.gz
        tabix -p vcf ~{prefix}.deletions.vcf.gz
    >>>

    output {
        File ins_vcf = "~{prefix}.insertions.vcf.gz"
        File ins_vcf_idx = "~{prefix}.insertions.vcf.gz.tbi"
        File del_vcf = "~{prefix}.deletions.vcf.gz"
        File del_vcf_idx = "~{prefix}.deletions.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 2 + 10,
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

task GenerateTRF {
    input {
        File vcf
        File vcf_idx
        String prefix
        String mode
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        if [[ ~{vcf} == *.gz ]]; then
            gunzip -c ~{vcf} > work_dir/input.vcf
            vcf_input="work_dir/input.vcf"
        else
            vcf_input="~{vcf}"
        fi

        if [[ "~{mode}" == "ins" ]]; then
            python3 /app/SVAN/scripts/ins2fasta.py "$vcf_input" work_dir
            trf work_dir/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.~{mode}_trf.out
        elif [[ "~{mode}" == "del" ]]; then
            python3 /app/SVAN/scripts/del2fasta.py "$vcf_input" work_dir
            trf work_dir/deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.~{mode}_trf.out
        else
            echo "Invalid mode provided to GenerateTRF task: ~{mode}"
            exit 1
        fi
    >>>

    output {
        File trf_output = "~{prefix}.~{mode}_trf.out"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 3 + 20,
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

task RunSvanAnnotation {
    input {
        File vcf
        File vcf_idx
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File mei_fa
        File mei_fa_amb
        File mei_fa_ann
        File mei_fa_bwt
        File mei_fa_pac
        File mei_fa_sa
        File mei_fa_mmi
        File ref_fa
        String prefix
        String mode
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fa, "GB") + size(mei_fa, "GB") + size(mei_fa_amb, "GB") + size(mei_fa_ann, "GB") + size(mei_fa_bwt, "GB") + size(mei_fa_pac, "GB") + size(mei_fa_sa, "GB") + size(mei_fa_mmi, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        if [[ ~{vcf} == *.gz ]]; then
            gunzip -c ~{vcf} > work_dir/input.vcf
            vcf_input="work_dir/input.vcf"
        else
            vcf_input="~{vcf}"
        fi
        
        svan_script_name=""
        if [[ "~{mode}" == "ins" ]]; then
            svan_script_name="SVAN-INS.py"
        elif [[ "~{mode}" == "del" ]]; then
            svan_script_name="SVAN-DEL.py"
        else
            echo "Invalid mode provided."
            exit 1
        fi

        python3 /app/SVAN/$svan_script_name \
            "$vcf_input" \
            ~{trf_output} \
            ~{vntr_bed} \
            ~{exons_bed} \
            ~{repeats_bed} \
            ~{mei_fa} \
            ~{ref_fa} \
            svan_annotated \
            -o work_dir

        bcftools sort work_dir/svan_annotated.vcf -Oz -o ~{prefix}.~{mode}.vcf.gz
        
        tabix -p vcf ~{prefix}.~{mode}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.~{mode}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.~{mode}.vcf.gz.tbi"
    }

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
