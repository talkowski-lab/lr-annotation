version 1.0

import "general/Structs.wdl"
import "general/Helpers.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_idx

        String prefix
        Int min_svlen

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

        String svan_docker
        Int? variants_per_shard

        RuntimeAttr? runtime_attr_separate
        RuntimeAttr? runtime_attr_generate_trf_ins
        RuntimeAttr? runtime_attr_generate_trf_del
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_merge
    }

    Int variants_per_shard_eff = select_first([variants_per_shard, 1000000000])

    call FilterBySvlen {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            min_svlen = min_svlen,
            prefix = "~{prefix}.filtered",
            docker = svan_docker
    }

    call Helpers.SplitVcfIntoShards {
        input:
            input_vcf = FilterBySvlen.filtered_vcf,
            input_vcf_idx = FilterBySvlen.filtered_vcf_idx,
            variants_per_shard = variants_per_shard_eff,
            prefix = "~{prefix}.split",
            docker_image = svan_docker
    }

    scatter (shard in zip(SplitVcfIntoShards.split_vcfs, SplitVcfIntoShards.split_vcf_indexes)) {
        String shard_prefix = basename(shard.left, ".vcf.gz")

        call SeparateInsertionsDeletions {
            input:
                vcf = shard.left,
                vcf_index = shard.right,
                prefix = shard_prefix,
                docker = svan_docker,
                runtime_attr_override = runtime_attr_separate
        }

        call GenerateTRF as GenerateTRFForInsertions {
            input:
                vcf = SeparateInsertionsDeletions.ins_vcf,
                vcf_index = SeparateInsertionsDeletions.ins_vcf_index,
                prefix = shard_prefix,
                mode = "ins",
                docker = svan_docker,
                runtime_attr_override = runtime_attr_generate_trf_ins
        }

        call GenerateTRF as GenerateTRFForDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_index = SeparateInsertionsDeletions.del_vcf_index,
                prefix = shard_prefix,
                mode = "del",
                docker = svan_docker,
                runtime_attr_override = runtime_attr_generate_trf_del
        }

        call RunSvanAnnotation as AnnotateInsertions {
            input:
                vcf = SeparateInsertionsDeletions.ins_vcf,
                vcf_index = SeparateInsertionsDeletions.ins_vcf_index,
                trf_output = GenerateTRFForInsertions.trf_output,
                vntr_bed = vntr_bed,
                exons_bed = exons_bed,
                repeats_bed = repeats_bed,
                mei_fa = mei_fa,
                mei_bwa_amb = mei_bwa_amb,
                mei_bwa_ann = mei_bwa_ann,
                mei_bwa_bwt = mei_bwa_bwt,
                mei_bwa_pac = mei_bwa_pac,
                mei_bwa_sa = mei_bwa_sa,
                mei_minimap_mmi = mei_minimap_mmi,
                ref_fa = ref_fa,
                prefix = shard_prefix,
                mode = "ins",
                docker = svan_docker,
                runtime_attr_override = runtime_attr_annotate_ins
        }

        call RunSvanAnnotation as AnnotateDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_index = SeparateInsertionsDeletions.del_vcf_index,
                trf_output = GenerateTRFForDeletions.trf_output,
                vntr_bed = vntr_bed,
                exons_bed = exons_bed,
                repeats_bed = repeats_bed,
                mei_fa = mei_fa,
                mei_bwa_amb = mei_bwa_amb,
                mei_bwa_ann = mei_bwa_ann,
                mei_bwa_bwt = mei_bwa_bwt,
                mei_bwa_pac = mei_bwa_pac,
                mei_bwa_sa = mei_bwa_sa,
                mei_minimap_mmi = mei_minimap_mmi,
                ref_fa = ref_fa,
                prefix = shard_prefix,
                mode = "del",
                docker = svan_docker,
                runtime_attr_override = runtime_attr_annotate_del
        }
    }
    
    call Helpers.ConcatVcfs as ConcatIns {
        input:
            vcfs = AnnotateInsertions.annotated_vcf,
            vcfs_idx = AnnotateInsertions.annotated_vcf_index,
            outfile_prefix = prefix + ".insertions",
            docker_image = svan_docker
    }

    call Helpers.ConcatVcfs as ConcatDel {
        input:
            vcfs = AnnotateDeletions.annotated_vcf,
            vcfs_idx = AnnotateDeletions.annotated_vcf_index,
            outfile_prefix = prefix + ".deletions",
            docker_image = svan_docker
    }

    call Helpers.ConcatVcfs as ConcatOther {
        input:
            vcfs = SeparateInsertionsDeletions.other_vcf,
            vcfs_idx = SeparateInsertionsDeletions.other_vcf_index,
            outfile_prefix = prefix + ".other",
            docker_image = svan_docker
    }

    call MergeAnnotatedVcfs {
        input:
            annotated_ins_vcf = ConcatIns.concat_vcf,
            annotated_ins_vcf_index = ConcatIns.concat_vcf_idx,
            annotated_del_vcf = ConcatDel.concat_vcf,
            annotated_del_vcf_index = ConcatDel.concat_vcf_idx,
            other_vcf = ConcatOther.concat_vcf,
            other_vcf_index = ConcatOther.concat_vcf_idx,
            prefix = prefix + "_annotated_filtered",
            docker = svan_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call MergeWithOriginal {
        input:
            original_vcf = vcf,
            original_vcf_idx= vcf_idx,
            annotated_vcf = MergeAnnotatedVcfs.merged_vcf,
            annotated_vcf_idx = MergeAnnotatedVcfs.merged_vcf_index,
            prefix = prefix,
            docker = svan_docker
    }

    output {
        File svan_annotated_vcf = MergeWithOriginal.final_vcf
        File svan_annotated_vcf_index = MergeWithOriginal.final_vcf_idx
    }
}

task FilterBySvlen {
    input {
        File vcf
        File vcf_idx
        String prefix
        Int min_svlen
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include 'abs(INFO/SVLEN) >= ~{min_svlen}' \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB") * 2) + 10,
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

task MergeWithOriginal {
    input {
        File original_vcf
        File original_vcf_idx
        File annotated_vcf
        File annotated_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -a ~{annotated_vcf} \
            -c "INFO" \
            -Oz -o "~{prefix}.svan_annotated.vcf.gz" \
            ~{original_vcf}

        tabix -p vcf "~{prefix}.svan_annotated.vcf.gz"
    >>>

    output {
        File final_vcf = "~{prefix}.svan_annotated.vcf.gz"
        File final_vcf_idx = "~{prefix}.svan_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(original_vcf, "GB") + size(annotated_vcf, "GB")) + 20,
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

task SeparateInsertionsDeletions {
    input {
        File vcf
        File vcf_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include 'INFO/SVTYPE="INS"' \
            -Oz -o ~{prefix}.insertions.vcf.gz
        
        bcftools view ~{vcf} \
            --include 'INFO/SVTYPE="DEL"' \
            -Oz -o ~{prefix}.deletions.vcf.gz
        
        bcftools view ~{vcf} \
            --exclude '(INFO/SVTYPE="INS" || INFO/SVTYPE="DEL")' \
            -Oz -o ~{prefix}.other_variants.vcf.gz
        
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
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB") * 2) + 10,
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

task GenerateTRF {
    input {
        File vcf
        File vcf_index
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
        mem_gb: 7.5,
        disk_gb: ceil(size(vcf, "GB") * 3) + 20,
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

task RunSvanAnnotation {
    input {
        File vcf
        File vcf_index
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File mei_fa
        File mei_bwa_amb
        File mei_bwa_ann
        File mei_bwa_bwt
        File mei_bwa_pac
        File mei_bwa_sa
        File mei_minimap_mmi
        File ref_fa
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
        
        svan_script_name=""
        if [[ "~{mode}" == "ins" ]]; then
            svan_script_name="SVAN-INS.py"
        elif [[ "~{mode}" == "del" ]]; then
            svan_script_name="SVAN-DEL.py"
        else
            echo "Invalid mode provided to RunSvanAnnotation task: ~{mode}"
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

        python3 <<'CODE'
import sys

def get_new_info_fields(svan_vcf_path):
    new_info_fields = {}
    with open(svan_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID='):
                parts = line.strip().split(',', 1)
                field_id = parts[0][len('##INFO=<ID='):]
                new_info_fields[field_id] = line.strip()
            elif line.startswith('#CHROM'):
                break
    return new_info_fields

def get_svan_annotations(svan_vcf_path):
    annotations = {}
    with open(svan_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, vcf_id, _, _, _, _, info = parts[:8]
            annotations[(chrom, pos, vcf_id)] = info
    return annotations

def merge_vcfs(original_vcf_path, svan_vcf_path, output_vcf_path):
    new_info_header_lines = get_new_info_fields(svan_vcf_path)
    svan_annotations = get_svan_annotations(svan_vcf_path)

    with open(original_vcf_path, 'r') as original_vcf, open(output_vcf_path, 'w') as out_vcf:
        original_info_fields = set()
        
        for line in original_vcf:
            if line.startswith('##INFO=<ID='):
                field_id = line.strip().split(',', 1)[0][len('##INFO=<ID='):]
                original_info_fields.add(field_id)
                out_vcf.write(line)
            elif line.startswith('##'):
                out_vcf.write(line)
            elif line.startswith('#CHROM'):
                for field_id, header_line in new_info_header_lines.items():
                    if field_id not in original_info_fields:
                        out_vcf.write(header_line + '\n')
                out_vcf.write(line)
                break
        
        for line in original_vcf:
            parts = line.strip().split('\t')
            chrom, pos, vcf_id, _, _, _, _, info = parts[:8]
            
            key = (chrom, pos, vcf_id)
            if key in svan_annotations:
                new_info = svan_annotations[key]
                if info == '.':
                    parts[7] = new_info
                else:
                    parts[7] = f"{info};{new_info}"
            
            out_vcf.write('\t'.join(parts) + '\n')

merge_vcfs('work_dir/input.vcf', 'work_dir/svan_annotated.vcf', 'work_dir/merged_annotated.vcf')
CODE

        bcftools sort work_dir/merged_annotated.vcf -Oz -o ~{prefix}.~{mode}_annotated.vcf.gz
        tabix -p vcf ~{prefix}.~{mode}_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.~{mode}_annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.~{mode}_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 15,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fa, "GB") + size(mei_fa, "GB") + size(mei_bwa_amb, "GB") + size(mei_bwa_ann, "GB") + size(mei_bwa_bwt, "GB") + size(mei_bwa_pac, "GB") + size(mei_bwa_sa, "GB") + size(mei_minimap_mmi, "GB")) + 25,
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

task MergeAnnotatedVcfs {
    input {
        File annotated_ins_vcf
        File annotated_ins_vcf_index
        File annotated_del_vcf
        File annotated_del_vcf_index
        File other_vcf
        File other_vcf_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        echo "~{annotated_ins_vcf}" > vcf_list.txt
        echo "~{annotated_del_vcf}" >> vcf_list.txt
        echo "~{other_vcf}" >> vcf_list.txt
        
        bcftools concat -a -f vcf_list.txt -Oz -o ~{prefix}.unsorted.vcf.gz
        
        bcftools sort ~{prefix}.unsorted.vcf.gz -Oz -o ~{prefix}.svan_annotated.vcf.gz
        
        tabix -p vcf ~{prefix}.svan_annotated.vcf.gz
        
        rm ~{prefix}.unsorted.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.svan_annotated.vcf.gz"
        File merged_vcf_index = "~{prefix}.svan_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 7.5,
        disk_gb: 50,
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