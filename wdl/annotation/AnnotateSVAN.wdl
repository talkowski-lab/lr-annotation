version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs

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

        String annotate_svan_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_separate
        RuntimeAttr? runtime_attr_generate_trf_ins
        RuntimeAttr? runtime_attr_generate_trf_del
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_concat
    }

    call FilterBySvlen {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            min_svlen = min_svlen,
            prefix = "~{prefix}.filtered",
            docker = annotate_svan_docker
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = FilterBySvlen.filtered_vcf,
                vcf_index = FilterBySvlen.filtered_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call SeparateInsertionsDeletions {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_index = SubsetVcfToContig.subset_vcf_index,
                prefix = "~{prefix}.~{contig}",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_separate
        }

        call GenerateTRF as GenerateTRFForInsertions {
            input:
                vcf = SeparateInsertionsDeletions.ins_vcf,
                vcf_index = SeparateInsertionsDeletions.ins_vcf_index,
                prefix = "~{prefix}.~{contig}",
                mode = "ins",
                docker = annotate_svan_docker,
                runtime_attr_override = runtime_attr_generate_trf_ins
        }

        call GenerateTRF as GenerateTRFForDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_index = SeparateInsertionsDeletions.del_vcf_index,
                prefix = "~{prefix}.~{contig}",
                mode = "del",
                docker = annotate_svan_docker,
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

        call RunSvanAnnotation as AnnotateDeletions {
            input:
                vcf = SeparateInsertionsDeletions.del_vcf,
                vcf_index = SeparateInsertionsDeletions.del_vcf_index,
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
    }
    
    call Helpers.ConcatTsvs as ConcatIns {
        input:
            tsvs = AnnotateInsertions.annotations_tsv,
            prefix = prefix + ".insertions",
            docker = annotate_svan_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call Helpers.ConcatTsvs as ConcatDel {
        input:
            tsvs = AnnotateDeletions.annotations_tsv,
            prefix = prefix + ".deletions",
            docker = annotate_svan_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call Helpers.ConcatTsvs as MergeFinal {
        input:
            tsvs = [ConcatIns.concatenated_tsv, ConcatDel.concatenated_tsv],
            prefix = prefix + ".svan_annotations",
            docker = annotate_svan_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call MergeHeaderFiles {
        input:
            header_files = flatten([AnnotateInsertions.header_file, AnnotateDeletions.header_file]),
            prefix = prefix,
            docker = annotate_svan_docker
    }

    output {
        File annotations_tsv_svan = MergeFinal.concatenated_tsv
        File annotations_header_svan = MergeHeaderFiles.merged_header
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
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 2 + 10,
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
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 2 + 10,
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
        disk_gb: ceil(size(vcf, "GB")) * 3 + 20,
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
        mem_gb: 15,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fa, "GB") + size(mei_fa, "GB") + size(mei_fa_amb, "GB") + size(mei_fa_ann, "GB") + size(mei_fa_bwt, "GB") + size(mei_fa_pac, "GB") + size(mei_fa_sa, "GB") + size(mei_fa_mmi, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

def get_new_info_field_ids(svan_vcf_path, original_vcf_path):
    original_info_fields = set()
    with open(original_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID='):
                field_id = line.strip().split(',', 1)[0][len('##INFO=<ID='):]
                original_info_fields.add(field_id)
            elif line.startswith('#CHROM'):
                break
    
    svan_info_fields = set()
    with open(svan_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID='):
                field_id = line.strip().split(',', 1)[0][len('##INFO=<ID='):]
                svan_info_fields.add(field_id)
            elif line.startswith('#CHROM'):
                break
    
    return sorted(svan_info_fields - original_info_fields)

def parse_info_field(info_string, field_ids):
    info_dict = {}
    if info_string == '.':
        return {field_id: '.' for field_id in field_ids}
    
    for item in info_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = 'TRUE'
    
    return {field_id: info_dict.get(field_id, '.') for field_id in field_ids}

def get_svan_annotations(svan_vcf_path, new_field_ids):
    annotations = {}
    with open(svan_vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, vcf_id, _, _, _, _, info = parts[:8]
            field_values = parse_info_field(info, new_field_ids)
            annotations[(chrom, pos, vcf_id)] = field_values
    return annotations

def create_tsv(original_vcf_path, svan_vcf_path, output_tsv_path, header_path):
    new_field_ids = get_new_info_field_ids(svan_vcf_path, original_vcf_path)
    svan_annotations = get_svan_annotations(svan_vcf_path, new_field_ids)
    
    # Write header file separately
    with open(header_path, 'w') as header_file:
        header_file.write('\t'.join(new_field_ids) + '\n')
    
    # Write TSV without header
    with open(output_tsv_path, 'w') as out_tsv:
        with open(original_vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom, pos, vcf_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                
                key = (chrom, pos, vcf_id)
                if key in svan_annotations:
                    field_values = svan_annotations[key]
                    row = [chrom, pos, ref, alt, vcf_id] + [field_values.get(field_id, '.') for field_id in new_field_ids]
                    out_tsv.write('\t'.join(row) + '\n')

create_tsv('work_dir/input.vcf', 'work_dir/svan_annotated.vcf', '~{prefix}.~{mode}_annotations.tsv', '~{prefix}.~{mode}_header.txt')
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.~{mode}_annotations.tsv"
        File header_file = "~{prefix}.~{mode}_header.txt"
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

task MergeHeaderFiles {
    input {
        Array[File] header_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Concatenate all header files and get unique fields while preserving order
        cat ~{sep=' ' header_files} | tr '\t' '\n' | awk '!seen[$0]++' | tr '\n' '\t' | sed 's/\t$/\n/' > ~{prefix}.svan_header.txt
    >>>

    output {
        File merged_header = "~{prefix}.svan_header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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