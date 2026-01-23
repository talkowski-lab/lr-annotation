version 1.0

import "../utils/Helpers.wdl" as Helpers
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

    call CheckSampleConsistency {
        input:
            snv_indel_vcf = snv_indel_vcf,
            sv_vcf = sv_vcf,
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

        call AnnotateVariantAttributes as AnnotateSnvIndel {
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

        call AnnotateVariantAttributes as AnnotateSv {
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

task CheckSampleConsistency {
    input {
        File snv_indel_vcf
        File sv_vcf
        Array[String] sample_ids
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' sample_ids} > requested_samples.txt

        bcftools query -l ~{snv_indel_vcf} | sort > snv_indel_samples.txt
        bcftools query -l ~{sv_vcf} | sort > sv_samples.txt

        missing_in_snv=""
        missing_in_sv=""
        
        while read sample; do
            if ! grep -q "^${sample}$" snv_indel_samples.txt; then
                missing_in_snv="${missing_in_snv} ${sample}"
            fi
            if ! grep -q "^${sample}$" sv_samples.txt; then
                missing_in_sv="${missing_in_sv} ${sample}"
            fi
        done < requested_samples.txt

        if [ -n "${missing_in_snv}" ]; then
            echo "ERROR: The following samples are missing from SNV/indel VCF:${missing_in_snv}"
            exit 1
        fi

        if [ -n "${missing_in_sv}" ]; then
            echo "ERROR: The following samples are missing from SV VCF:${missing_in_sv}"
            exit 1
        fi

        echo "All requested samples are present in both VCFs"
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([snv_indel_vcf, sv_vcf], "GB")) + 10,
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
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task AnnotateVariantAttributes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=VARLEN'; then
            echo '##INFO=<ID=VARLEN,Number=1,Type=Integer,Description="Variant length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=VARTYPE'; then
            echo '##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVLEN'; then
            echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVTYPE'; then
            echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">' >> new_headers.txt
        fi

        bcftools annotate \
            -h new_headers.txt \
            ~{vcf} \
            -Oz -o temp.vcf.gz
        tabix -p vcf temp.vcf.gz

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VARLEN\t%INFO/VARTYPE\t%INFO/SVLEN\t%INFO/SVTYPE\n' temp.vcf.gz | \
        awk -F'\t' '{
            split($4, alleles, ",")
            ref_len = length($3)
            best_type = "SNV"
            max_diff = 0
            
            alt_len = length($4)
            diff = (alt_len > ref_len) ? (alt_len - ref_len) : (ref_len - alt_len)
            max_diff = diff
            if (alt_len > ref_len) {
                best_type = "INS"
            } else if (alt_len < ref_len) {
                best_type = "DEL"
            } else {
                best_type = "SNV"
            }

            var_len = max_diff
            var_type = best_type

            if (var_type == "INS" || var_type == "DEL") {
                sv_len = var_len
                sv_type = var_type
            } else {
                var_len = ($5 != ".") ? $5 : "."
                var_type = ($6 != ".") ? $6 : "."
                sv_len = ($7 != ".") ? $7 : "."
                sv_type = ($8 != ".") ? $8 : "."
            }

            print $1"\t"$2"\t"$3"\t"$4"\t"var_len"\t"var_type"\t"sv_len"\t"sv_type
        }' \
            | bgzip -c > annot.txt.gz
        
        tabix -s1 -b2 -e2 annot.txt.gz

        bcftools annotate -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/VARLEN,INFO/VARTYPE,INFO/SVLEN,INFO/SVTYPE \
            temp.vcf.gz \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 5,
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

