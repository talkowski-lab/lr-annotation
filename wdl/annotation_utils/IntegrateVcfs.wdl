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
        Int? records_per_shard
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
        RuntimeAttr? runtime_attr_shard_snv_indel
        RuntimeAttr? runtime_attr_split_snv_indel
        RuntimeAttr? runtime_attr_subset_samples_snv_indel
        RuntimeAttr? runtime_attr_annotate_attributes_snv_indel
        RuntimeAttr? runtime_attr_add_info_snv_indel
        RuntimeAttr? runtime_attr_add_filter_snv_indel
        RuntimeAttr? runtime_attr_concat_snv_indel_shards

        RuntimeAttr? runtime_attr_subset_contig_sv
        RuntimeAttr? runtime_attr_shard_sv
        RuntimeAttr? runtime_attr_split_sv
        RuntimeAttr? runtime_attr_subset_samples_sv
        RuntimeAttr? runtime_attr_annotate_attributes_sv
        RuntimeAttr? runtime_attr_add_info_sv
        RuntimeAttr? runtime_attr_add_filter_sv
        RuntimeAttr? runtime_attr_concat_sv_shards
        
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_rename_and_filter
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(swap_samples_snv_indel)) {
        call Helpers.SwapSampleIds as SwapSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                sample_swap_list = select_first([swap_samples_snv_indel]),
                prefix = "~{prefix}.snv_indel.swapped",
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
                prefix = "~{prefix}.sv.swapped",
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

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords as ShardSnvIndel {
                input:
                    vcf = SubsetContigSnvIndel.subset_vcf,
                    vcf_idx = SubsetContigSnvIndel.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.snv_indel",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard_snv_indel
            }
        }

        Array[File] snv_indel_vcfs_to_process = select_first([ShardSnvIndel.shards, [SubsetContigSnvIndel.subset_vcf]])
        Array[File] snv_indel_vcf_idxs_to_process = select_first([ShardSnvIndel.shard_idxs, [SubsetContigSnvIndel.subset_vcf_idx]])

        scatter (i in range(length(snv_indel_vcfs_to_process))) {
            call Helpers.SplitMultiallelics as SplitSnvIndel {
                input:
                    vcf = snv_indel_vcfs_to_process[i],
                    vcf_idx = snv_indel_vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.split",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_snv_indel
            }

            call Helpers.SubsetVcfToSamples as SubsetSamplesSnvIndel {
                input:
                    vcf = SplitSnvIndel.split_vcf,
                    vcf_idx = SplitSnvIndel.split_vcf_idx,
                    samples = sample_ids,
                    prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_samples_snv_indel
            }

            call Helpers.AnnotateVariantAttributes as AnnotateSnvIndel {
                input:
                    vcf = SubsetSamplesSnvIndel.subset_vcf,
                    vcf_idx = SubsetSamplesSnvIndel.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.annotated",
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
                    prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.add_info",
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
                    prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.add_filter",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_add_filter_snv_indel
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatSnvIndelShards {
                input:
                    vcfs = AddFilterSnvIndel.flagged_vcf,
                    vcf_idxs = AddFilterSnvIndel.flagged_vcf_idx,
                    prefix = "~{prefix}.~{contig}.snv_indel.concatenated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_snv_indel_shards
            }
        }

        File final_snv_indel_vcf_for_contig = select_first([ConcatSnvIndelShards.concat_vcf, AddFilterSnvIndel.flagged_vcf[0]])
        File final_snv_indel_vcf_idx_for_contig = select_first([ConcatSnvIndelShards.concat_vcf_idx, AddFilterSnvIndel.flagged_vcf_idx[0]])

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

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords as ShardSv {
                input:
                    vcf = SubsetContigSv.subset_vcf,
                    vcf_idx = SubsetContigSv.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.sv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard_sv
            }
        }

        Array[File] sv_vcfs_to_process = select_first([ShardSv.shards, [SubsetContigSv.subset_vcf]])
        Array[File] sv_vcf_idxs_to_process = select_first([ShardSv.shard_idxs, [SubsetContigSv.subset_vcf_idx]])

        scatter (i in range(length(sv_vcfs_to_process))) {
            call Helpers.SplitMultiallelics as SplitSv {
                input:
                    vcf = sv_vcfs_to_process[i],
                    vcf_idx = sv_vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.sv.shard_~{i}.split",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_sv
            }

            call Helpers.SubsetVcfToSamples as SubsetSamplesSv {
                input:
                    vcf = SplitSv.split_vcf,
                    vcf_idx = SplitSv.split_vcf_idx,
                    samples = sample_ids,
                    prefix = "~{prefix}.~{contig}.sv.shard_~{i}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_samples_sv
            }

            call Helpers.AnnotateVariantAttributes as AnnotateSv {
                input:
                    vcf = SubsetSamplesSv.subset_vcf,
                    vcf_idx = SubsetSamplesSv.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.sv.shard_~{i}.annotated",
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
                    prefix = "~{prefix}.~{contig}.sv.shard_~{i}.add_info",
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
                    prefix = "~{prefix}.~{contig}.sv.shard_~{i}.add_filter",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_add_filter_sv
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatSvShards {
                input:
                    vcfs = AddFilterSv.flagged_vcf,
                    vcf_idxs = AddFilterSv.flagged_vcf_idx,
                    prefix = "~{prefix}.~{contig}.sv.concatenated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_sv_shards
            }
        }

        File final_sv_vcf_for_contig = select_first([ConcatSvShards.concat_vcf, AddFilterSv.flagged_vcf[0]])
        File final_sv_vcf_idx_for_contig = select_first([ConcatSvShards.concat_vcf_idx, AddFilterSv.flagged_vcf_idx[0]])

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [final_snv_indel_vcf_for_contig, final_sv_vcf_for_contig],
                vcf_idxs = [final_snv_indel_vcf_idx_for_contig, final_sv_vcf_idx_for_contig],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call RenameAndFilterVariants {
            input:
                vcf = MergeContigVcfs.concat_vcf,
                vcf_idx = MergeContigVcfs.concat_vcf_idx,
                snv_indel_vcf_size_flag = snv_indel_vcf_size_flag,
                snv_indel_vcf_source_tag = snv_indel_vcf_source_tag,
                sv_vcf_source_tag = sv_vcf_source_tag,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_and_filter
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = RenameAndFilterVariants.filtered_vcf,
            vcf_idxs = RenameAndFilterVariants.filtered_vcf_idx,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task RenameAndFilterVariants {
    input {
        File vcf
        File vcf_idx
        String snv_indel_vcf_size_flag
        String snv_indel_vcf_source_tag
        String sv_vcf_source_tag
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile
from collections import defaultdict

snv_indel_source = "~{snv_indel_vcf_source_tag}"
sv_source = "~{sv_vcf_source_tag}"
size_flag = "~{snv_indel_vcf_size_flag}"

# First pass: collect SV variants for redundancy check
vcf_in = VariantFile("~{vcf}")
sv_variants = set()
for record in vcf_in:
    if record.info.get("SOURCE") == sv_source:
        sv_variants.add((record.chrom, record.pos, record.ref, tuple(record.alts)))
vcf_in.close()

# Second pass: collect ID counts for non-redundant variants
vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    # Skip redundant variants
    if (record.info.get("SOURCE") == snv_indel_source and 
        size_flag in record.filter and 
        (record.chrom, record.pos, record.ref, tuple(record.alts)) in sv_variants):
        continue
    
    # Rename variant IDs
    a_type = record.info.get('allele_type').upper()
    a_len = abs(int(record.info.get('allele_length')))
    if a_type == "SNV":
        new_id = f"{record.chrom}-{record.pos}-{a_type}"
    else:
        new_id = f"{record.chrom}-{record.pos}-{a_type}-{a_len}"
    
    id_counts[new_id] += 1
vcf_in.close()

# Third pass: rename IDs and filter redundant variants
vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    # Skip redundant variants
    if (record.info.get("SOURCE") == snv_indel_source and 
        size_flag in record.filter and 
        (record.chrom, record.pos, record.ref, tuple(record.alts)) in sv_variants):
        continue
    
    # Generate new ID
    a_type = record.info.get('allele_type').upper()
    a_len = abs(int(record.info.get('allele_length')))
    if a_type == "SNV":
        new_id = f"{record.chrom}-{record.pos}-{a_type}"
    else:
        new_id = f"{record.chrom}-{record.pos}-{a_type}-{a_len}"
    
    # Handle duplicate IDs
    if id_counts[new_id] > 1:
        id_seen[new_id] += 1
        record.id = f"{new_id}_{id_seen[new_id]}"
    else:
        record.id = new_id
    
    vcf_out.write(record)
vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
