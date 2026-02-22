version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        
        Int? records_per_shard
        
        File vntr_bed
        File exons_bed
        File repeats_bed
        File mei_fa
        Array[File] mei_fa_indices
        File ref_fa
        Array[File] ref_fa_indices

        String utils_docker
        String svan_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_reset_filters
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_subset_ins
        RuntimeAttr? runtime_attr_subset_del
        RuntimeAttr? runtime_attr_generate_trf_ins
        RuntimeAttr? runtime_attr_generate_trf_del
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_extract_ins
        RuntimeAttr? runtime_attr_extract_del
        RuntimeAttr? runtime_attr_concat_ins
        RuntimeAttr? runtime_attr_concat_del
        RuntimeAttr? runtime_attr_concat_final
    }

    scatter (contig in contigs) {
        # Preprocessing
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                extra_args = "-G",
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ResetVcfFilters {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.reset_filters",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_reset_filters
        }

        # Insertions
        call Helpers.SubsetVcfByArgs as SubsetIns {
            input:
                vcf = ResetVcfFilters.reset_vcf,
                vcf_idx = ResetVcfFilters.reset_vcf_idx,
                include_args = 'INFO/allele_type="ins"',
                prefix = "~{prefix}.~{contig}.ins_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_ins
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords as ShardIns {
                input:
                    vcf = SubsetIns.subset_vcf,
                    vcf_idx = SubsetIns.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.ins",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] ins_vcfs_to_process = select_first([ShardIns.shards, [SubsetIns.subset_vcf]])
        Array[File] ins_vcf_idxs_to_process = select_first([ShardIns.shard_idxs, [SubsetIns.subset_vcf_idx]])

        scatter (ins_shard_idx in range(length(ins_vcfs_to_process))) {
            call GenerateTRF as GenerateTRFIns {
                input:
                    vcf = ins_vcfs_to_process[ins_shard_idx],
                    vcf_idx = ins_vcf_idxs_to_process[ins_shard_idx],
                    prefix = "~{prefix}.~{contig}.ins_shard_~{ins_shard_idx}.trf",
                    mode = "ins",
                    docker = svan_docker,
                    runtime_attr_override = runtime_attr_generate_trf_ins
            }

            call RunSvanAnnotate as SvanAnnotateIns {
                input:
                    vcf = ins_vcfs_to_process[ins_shard_idx],
                    vcf_idx = ins_vcf_idxs_to_process[ins_shard_idx],
                    trf_output = GenerateTRFIns.trf_output,
                    vntr_bed = vntr_bed,
                    exons_bed = exons_bed,
                    repeats_bed = repeats_bed,
                    mei_fa = mei_fa,
                    mei_fa_indices = mei_fa_indices,
                    ref_fa = ref_fa,
                    ref_fa_indices = ref_fa_indices,
                    prefix = "~{prefix}.~{contig}.ins_shard_~{ins_shard_idx}.svan",
                    mode = "ins",
                    docker = svan_docker,
                    runtime_attr_override = runtime_attr_annotate_ins
            }

            call Helpers.ExtractVcfAnnotations as ExtractIns {
                input:
                    vcf = SvanAnnotateIns.annotated_vcf,
                    vcf_idx = SvanAnnotateIns.annotated_vcf_idx,
                    original_vcf = ins_vcfs_to_process[ins_shard_idx],
                    original_vcf_idx = ins_vcf_idxs_to_process[ins_shard_idx],
                    add_header_row = true,
                    prefix = "~{prefix}.~{contig}.ins_shard_~{ins_shard_idx}.annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_ins
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatInsShards {
                input:
                    tsvs = ExtractIns.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.ins_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_ins
            }
        }

        File final_ins_annotations = select_first([ConcatInsShards.concatenated_tsv, ExtractIns.annotations_tsv[0]])

        # Deletions
        call Helpers.SubsetVcfByArgs as SubsetDel {
            input:
                vcf = ResetVcfFilters.reset_vcf,
                vcf_idx = ResetVcfFilters.reset_vcf_idx,
                include_args = 'INFO/allele_type="del"',
                prefix = "~{prefix}.~{contig}.del_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_del
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords as ShardDel {
                input:
                    vcf = SubsetDel.subset_vcf,
                    vcf_idx = SubsetDel.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.del",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] del_vcfs_to_process = select_first([ShardDel.shards, [SubsetDel.subset_vcf]])
        Array[File] del_vcf_idxs_to_process = select_first([ShardDel.shard_idxs, [SubsetDel.subset_vcf_idx]])

        scatter (del_shard_idx in range(length(del_vcfs_to_process))) {
            call GenerateTRF as GenerateTRFDel {
                input:
                    vcf = del_vcfs_to_process[del_shard_idx],
                    vcf_idx = del_vcf_idxs_to_process[del_shard_idx],
                    prefix = "~{prefix}.~{contig}.del_shard_~{del_shard_idx}.trf",
                    mode = "del",
                    docker = svan_docker,
                    runtime_attr_override = runtime_attr_generate_trf_del
            }

            call RunSvanAnnotate as SvanAnnotateDel {
                input:
                    vcf = del_vcfs_to_process[del_shard_idx],
                    vcf_idx = del_vcf_idxs_to_process[del_shard_idx],
                    trf_output = GenerateTRFDel.trf_output,
                    vntr_bed = vntr_bed,
                    exons_bed = exons_bed,
                    repeats_bed = repeats_bed,
                    mei_fa = mei_fa,
                    mei_fa_indices = mei_fa_indices,
                    ref_fa = ref_fa,
                    ref_fa_indices = ref_fa_indices,
                    prefix = "~{prefix}.~{contig}.del_shard_~{del_shard_idx}.svan",
                    mode = "del",
                    docker = svan_docker,
                    runtime_attr_override = runtime_attr_annotate_del
            }

            call Helpers.ExtractVcfAnnotations as ExtractDel {
                input:
                    vcf = SvanAnnotateDel.annotated_vcf,
                    vcf_idx = SvanAnnotateDel.annotated_vcf_idx,
                    original_vcf = del_vcfs_to_process[del_shard_idx],
                    original_vcf_idx = del_vcf_idxs_to_process[del_shard_idx],
                    add_header_row = true,
                    prefix = "~{prefix}.~{contig}.del_shard_~{del_shard_idx}.annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_del
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatDelShards {
                input:
                    tsvs = ExtractDel.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.del_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_del
            }
        }

        File final_del_annotations = select_first([ConcatDelShards.concatenated_tsv, ExtractDel.annotations_tsv[0]])
    }
    
    # Postprocessing
    call Helpers.ConcatAlignedTsvs {
        input:
            tsvs = flatten([final_ins_annotations, final_del_annotations]),
            prefix = prefix + ".svan_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_final
    }

    output {
        File annotations_tsv_svan = ConcatAlignedTsvs.merged_tsv
        File annotations_header_svan = ConcatAlignedTsvs.merged_header
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

        gunzip -c ~{vcf} > tmp.vcf

        mkdir -p work_dir

        if [[ "~{mode}" == "ins" ]]; then
            python3 /app/SVAN/scripts/ins2fasta.py tmp.vcf work_dir
            trf work_dir/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.out
        elif [[ "~{mode}" == "del" ]]; then
            python3 /app/SVAN/scripts/del2fasta.py tmp.vcf work_dir
            trf work_dir/deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.out
        fi
    >>>

    output {
        File trf_output = "~{prefix}.out"
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

task RunSvanAnnotate {
    input {
        File vcf
        File vcf_idx
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File mei_fa
        Array[File] mei_fa_indices
        File ref_fa
        Array[File] ref_fa_indices
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

        bcftools sort work_dir/svan_annotated.vcf -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(mei_fa, "GB") + size(mei_fa_indices, "GB") + size(ref_fa, "GB")  + size(ref_fa_indices, "GB")) + 20,
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
