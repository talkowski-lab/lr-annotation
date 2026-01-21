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
        Array[String] sample_ids
        Int max_size_snv_indel = 49
        Int min_size_sv = 50
        String prefix
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
        call Helpers.SubsetVcfToSampleList as SubsetSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_index = snv_indel_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv_indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call AnnotateSvlenSvtype as AnnotateSnvIndel {
            input:
                vcf = SubsetSnvIndel.subset_vcf,
                vcf_idx = SubsetSnvIndel.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.snv_indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.SubsetVcfBySize as FilterSnvIndel {
            input:
                vcf = AnnotateSnvIndel.annotated_vcf,
                vcf_index = AnnotateSnvIndel.annotated_vcf_idx,
                locus = contig,
                max_size = max_size_snv_indel,
                prefix = "~{prefix}.~{contig}.snv_indel",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = sv_vcf,
                vcf_index = sv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call AnnotateSvlenSvtype as AnnotateSv {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_index,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.SubsetVcfBySize as FilterSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_index = AnnotateSv.annotated_vcf_idx,
                locus = contig,
                min_size = min_size_sv,
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [FilterSnvIndel.subset_vcf, FilterSv.subset_vcf],
                vcfs_idx = [FilterSnvIndel.subset_vcf_idx, FilterSv.subset_vcf_idx],
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

task AnnotateSvlenSvtype {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Add missing headers so we can always query both fields
        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVLEN'; then
            echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVTYPE'; then
            echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">' >> new_headers.txt
        fi

        if [ -s new_headers.txt ]; then
            bcftools annotate -h new_headers.txt ~{vcf} -Oz -o temp.vcf.gz
            tabix -p vcf temp.vcf.gz
        else
            ln -s ~{vcf} temp.vcf.gz
            ln -s ~{vcf_idx} temp.vcf.gz.tbi
        fi

        # Query and annotate with unified awk script
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/SVTYPE\n' temp.vcf.gz | \
        awk '{
            ref_len = length($3)
            alt_len = length($4)
            
            if (ref_len == alt_len) {
                calc_svtype = (ref_len == 1) ? "SNV" : "INS"
                calc_svlen = (ref_len == 1) ? 1 : ref_len
            } else if (ref_len > alt_len) {
                calc_svtype = "DEL"
                calc_svlen = ref_len - alt_len
            } else {
                calc_svtype = "INS"
                calc_svlen = alt_len - ref_len
            }
            
            final_svtype = ($6 != ".") ? $6 : calc_svtype
            final_svlen = (ref_len == 1 && alt_len == 1) ? 1 : (($5 != ".") ? $5 : calc_svlen)
            
            print $1"\t"$2"\t"$3"\t"$4"\t"final_svtype"\t"final_svlen
        }' | bgzip -c > annot.txt.gz
        
        tabix -s1 -b2 -e2 annot.txt.gz

        bcftools annotate -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/SVTYPE,INFO/SVLEN \
            temp.vcf.gz -Oz -o ~{prefix}.vcf.gz

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

