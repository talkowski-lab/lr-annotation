version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateVcfs {
    input {
        File snv_indel_vcf
        File snv_indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        Int min_sv_length = 50
        String snv_indel_vcf_source_tag
        String snv_indel_vcf_size_flag
        String snv_indel_vcf_size_flag_description
        String sv_vcf_source_tag
        String sv_vcf_size_flag
        String sv_vcf_size_flag_description
        File? swap_samples_snv_indel
        File? swap_samples_sv
        
        String utils_docker
        
        RuntimeAttr? runtime_attr_swap_samples_snv_indel
        RuntimeAttr? runtime_attr_swap_samples_sv
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_contig_snv_indel
        RuntimeAttr? runtime_attr_split_snv_indel
        RuntimeAttr? runtime_attr_subset_samples_snv_indel
        RuntimeAttr? runtime_attr_annotate_attributes_snv_indel
        RuntimeAttr? runtime_attr_add_info_snv_indel
        RuntimeAttr? runtime_attr_add_filter_snv_indel
        RuntimeAttr? runtime_attr_rename_snv_indel
        RuntimeAttr? runtime_attr_subset_contig_sv
        RuntimeAttr? runtime_attr_split_sv
        RuntimeAttr? runtime_attr_subset_samples_sv
        RuntimeAttr? runtime_attr_annotate_attributes_sv
        RuntimeAttr? runtime_attr_add_info_sv
        RuntimeAttr? runtime_attr_add_filter_sv
        RuntimeAttr? runtime_attr_rename_sv
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(swap_samples_snv_indel)) {
        call Helpers.SwapSampleIds as SwapSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                sample_swap_list = select_first([swap_samples_snv_indel]),
                prefix = prefix + ".snv_indel.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples_snv_indel
        }
    }

    if (defined(swap_samples_sv)) {
        call Helpers.SwapSampleIds as SwapSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                sample_swap_list = select_first([swap_samples_sv]),
                prefix = prefix + ".sv.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples_sv
        }
    }
    
    File final_snv_indel_vcf = select_first([SwapSnvIndel.swapped_vcf, snv_indel_vcf])
    File final_snv_indel_vcf_idx = select_first([SwapSnvIndel.swapped_vcf_idx, snv_indel_vcf_idx])
    File final_sv_vcf = select_first([SwapSv.swapped_vcf, sv_vcf])
    File final_sv_vcf_idx = select_first([SwapSv.swapped_vcf_idx, sv_vcf_idx])

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = [final_snv_indel_vcf, final_sv_vcf],
            vcf_idxs = [final_snv_indel_vcf, final_sv_vcf],
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        # SNV Indel Processing
        call Helpers.SubsetVcfToContig as SubsetContigSnvIndel {
            input:
                vcf = final_snv_indel_vcf,
                vcf_idx = final_snv_indel_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv_indel",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_snv_indel
        }

        call Helpers.SplitMultiallelics as SplitSnvIndel {
            input:
                vcf = SubsetContigSnvIndel.subset_vcf,
                vcf_idx = SubsetContigSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_split_snv_indel
        }

        call Helpers.SubsetVcfToSampleList as SubsetSamplesSnvIndel {
            input:
                vcf = SplitSnvIndel.split_vcf,
                vcf_idx = SplitSnvIndel.split_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.snv_indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_snv_indel
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSnvIndel {
            input:
                vcf = SubsetSamplesSnvIndel.subset_vcf,
                vcf_idx = SubsetSamplesSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_attributes_snv_indel
        }

        call Helpers.AddInfo as AddInfoSnvIndel {
            input:
                vcf = AnnotateSnvIndel.annotated_vcf,
                vcf_idx = AnnotateSnvIndel.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = snv_indel_vcf_source_tag,
                tag_description = "Source of variant call",
                prefix = "~{prefix}.~{contig}.snv_indel.add_info",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info_snv_indel
        }

        call Helpers.AddFilter as AddFilterSnvIndel {
            input:
                vcf = AddInfoSnvIndel.annotated_vcf,
                vcf_idx = AddInfoSnvIndel.annotated_vcf_idx,
                filter_name = snv_indel_vcf_size_flag,
                filter_description = snv_indel_vcf_size_flag_description,
                filter_expression = "abs(INFO/allele_length) >= ~{min_sv_length}",
                prefix = "~{prefix}.~{contig}.snv_indel.add_filter",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_filter_snv_indel
        }

        call RenameSnvIndelIds {
            input:
                vcf = AddFilterSnvIndel.flagged_vcf,
                vcf_idx = AddFilterSnvIndel.flagged_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.renamed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_snv_indel
        }

        # SV Processing
        call Helpers.SubsetVcfToContig as SubsetContigSv {
            input:
                vcf = final_sv_vcf,
                vcf_idx = final_sv_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_sv
        }

        call Helpers.SplitMultiallelics as SplitSv {
            input:
                vcf = SubsetContigSv.subset_vcf,
                vcf_idx = SubsetContigSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_split_sv
        }

        call Helpers.SubsetVcfToSampleList as SubsetSamplesSv {
            input:
                vcf = SplitSv.split_vcf,
                vcf_idx = SplitSv.split_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_sv
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSv {
            input:
                vcf = SubsetSamplesSv.subset_vcf,
                vcf_idx = SubsetSamplesSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_attributes_sv
        }

        call Helpers.AddInfo as AddInfoSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_idx = AnnotateSv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = sv_vcf_source_tag,
                tag_description = "Source of variant call",
                prefix = "~{prefix}.~{contig}.sv.add_info",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info_sv
        }

        call Helpers.AddFilter as AddFilterSv {
            input:
                vcf = AddInfoSv.annotated_vcf,
                vcf_idx = AddInfoSv.annotated_vcf_idx,
                filter_name = sv_vcf_size_flag,
                filter_description = sv_vcf_size_flag_description,
                filter_expression = "abs(INFO/allele_length) < ~{min_sv_length}",
                prefix = "~{prefix}.~{contig}.sv.add_filter",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_filter_sv
        }

        call RenameSvIds {
            input:
                vcf = AddFilterSv.flagged_vcf,
                vcf_idx = AddFilterSv.flagged_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.renamed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_sv
        }

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [RenameSnvIndelIds.renamed_vcf, RenameSvIds.renamed_vcf],
                vcf_idxs = [RenameSnvIndelIds.renamed_vcf_idx, RenameSvIds.renamed_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcf_idxs = MergeContigVcfs.concat_vcf_idx,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task RenameSnvIndelIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile

vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)

for record in vcf_in:
    allele_type = record.info.get('allele_type').upper()
    allele_length = abs(int(record.info.get('allele_length')))
    record.id = f"{record.chrom}-{record.pos}-{allele_type}-{allele_length}"
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 8 * ceil(size(vcf, "GB")) + 5,
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

task RenameSvIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            --set-id '%CHROM-%POS-%REF-%ALT' \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 8 * ceil(size(vcf, "GB")) + 5,
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
