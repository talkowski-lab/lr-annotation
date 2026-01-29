version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateHPRC {
    input {
        File snv_indel_vcf
        File snv_indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        Int min_size_sv = 50
        String snv_indel_vcf_source_tag
        String snv_indel_vcf_source_tag_description
        String snv_indel_vcf_size_flag
        String snv_indel_vcf_size_flag_description
        String sv_vcf_source_tag
        String sv_vcf_source_tag_description
        String sv_vcf_size_flag
        String sv_vcf_size_flag_description
        
        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_filter_snv_indel
        RuntimeAttr? runtime_attr_filter_sv
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_annotate_svlen_svtype
    }

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = [snv_indel_vcf, sv_vcf],
            vcfs_idx = [snv_indel_vcf_idx, sv_vcf_idx],
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        # SNV Indel VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv_indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call SplitMultiallelics {
            input:
                vcf = SubsetSnvIndel.subset_vcf,
                vcf_idx = SubsetSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSnvIndel {
            input:
                vcf = SplitMultiallelics.split_vcf,
                vcf_idx = SplitMultiallelics.split_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.AddInfo as AddInfoSnvIndel {
            input:
                vcf = AnnotateSnvIndel.annotated_vcf,
                vcf_idx = AnnotateSnvIndel.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = snv_indel_vcf_source_tag,
                tag_description = snv_indel_vcf_source_tag_description,
                prefix = "~{prefix}.~{contig}.snv_indel.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.AddFilter as AddFilterSnvIndel {
            input:
                vcf = AddInfoSnvIndel.annotated_vcf,
                vcf_idx = AddInfoSnvIndel.annotated_vcf_idx,
                filter_name = snv_indel_vcf_size_flag,
                filter_description = snv_indel_vcf_size_flag_description,
                filter_expression = "abs(INFO/VARLEN) >= ~{min_size_sv}",
                prefix = "~{prefix}.~{contig}.snv_indel.flagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        # SV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call SplitMultiallelics as SplitSv {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSv {
            input:
                vcf = SplitSv.split_vcf,
                vcf_idx = SplitSv.split_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.AddInfo as AddInfoSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_idx = AnnotateSv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = sv_vcf_source_tag,
                tag_description = sv_vcf_source_tag_description,
                prefix = "~{prefix}.~{contig}.sv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.AddFilter as AddFilterSv {
            input:
                vcf = AddInfoSv.annotated_vcf,
                vcf_idx = AddInfoSv.annotated_vcf_idx,
                filter_name = sv_vcf_size_flag,
                filter_description = sv_vcf_size_flag_description,
                filter_expression = "abs(INFO/VARLEN) < ~{min_size_sv}",
                prefix = "~{prefix}.~{contig}.sv.flagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [AddFilterSnvIndel.flagged_vcf, AddFilterSv.flagged_vcf],
                vcfs_idx = [AddFilterSnvIndel.flagged_vcf_idx, AddFilterSv.flagged_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcfs_idx = MergeContigVcfs.concat_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task SplitMultiallelics {
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
import pysam
import sys

vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("temp.vcf", "w", header=vcf_in.header)

for record in vcf_in:
    if len(record.alts) <= 1:
        vcf_out.write(record)
        continue
    
    ids = record.id.split(';')
    for i, alt_seq in enumerate(record.alts):
        parts = ids[i].split('_')        
        new_rec = record.copy()
        new_rec.chrom = parts[0]
        new_rec.pos = int(parts[1])
        new_rec.ref = parts[2]
        new_rec.alts = (parts[3],)
        new_rec.id = ids[i]
        target_allele_idx = i + 1
        
        for sample in record.samples:
            old_gt = record.samples[sample]['GT']
            new_gt = []
            for allele in old_gt:
                if allele is None:
                    new_gt.append(None)
                elif allele == target_allele_idx:
                    new_gt.append(1)
                else:
                    new_gt.append(0)
            new_rec.samples[sample]['GT'] = tuple(new_gt)
        
        vcf_out.write(new_rec)

vcf_in.close()
vcf_out.close()
CODE

        bcftools sort temp.vcf -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File split_vcf = "~{prefix}.vcf.gz"
        File split_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(vcf, "GB")) + 5,
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
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
