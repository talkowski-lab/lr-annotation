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
        Int max_size_snv_indel = 50
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
        call SubsetAndFilterVcf as SubsetSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                contig = contig,
                sample_ids = sample_ids,
                size_threshold = max_size_snv_indel,
                use_max_size = true,
                prefix = "~{prefix}.~{contig}.snv_indel",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call SubsetAndFilterVcf as SubsetSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                contig = contig,
                sample_ids = sample_ids,
                size_threshold = min_size_sv,
                use_max_size = false,
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [SubsetSnvIndel.filtered_vcf, SubsetSv.filtered_vcf],
                vcfs_idx = [SubsetSnvIndel.filtered_vcf_idx, SubsetSv.filtered_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call AnnotateSvlenSvtype {
            input:
                vcf = MergeContigVcfs.concat_vcf,
                vcf_idx = MergeContigVcfs.concat_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateSvlenSvtype.annotated_vcf,
            vcfs_idx = AnnotateSvlenSvtype.annotated_vcf_idx,
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

task SubsetAndFilterVcf {
    input {
        File vcf
        File vcf_idx
        String contig
        Array[String] sample_ids
        Int size_threshold
        Boolean use_max_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' sample_ids} > samples.txt

        bcftools view -S samples.txt -c 1 -r ~{contig} ~{vcf} -Oz -o subset.vcf.gz
        tabix -p vcf subset.vcf.gz

        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' subset.vcf.gz | \
        awk -v threshold=~{size_threshold} -v use_max=~{use_max_size} '{
            ref_len = length($4)
            alt_len = length($5)
            size = (ref_len > alt_len) ? (ref_len - alt_len) : (alt_len - ref_len)
            if ((use_max == "true" && size <= threshold) || (use_max == "false" && size >= threshold)) {
                print $3
            }
        }' > keep_ids.txt

        if [ -s keep_ids.txt ]; then
            bcftools view -i 'ID=@keep_ids.txt' subset.vcf.gz -Oz -o ~{prefix}.vcf.gz
        else
            bcftools view -h subset.vcf.gz | bgzip -c > ~{prefix}.vcf.gz
        fi

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

        bcftools annotate -x INFO/SVLEN,INFO/SVTYPE ~{vcf} | \
        bcftools +fill-tags -Oz -o annotated.vcf.gz -- -t 'SVLEN:1=int(strlen(ALT[0])-strlen(REF))'
        tabix -p vcf annotated.vcf.gz

        touch headers.txt
        if ! bcftools view -h annotated.vcf.gz | grep -q '##INFO=<ID=SVLEN'; then
            echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">' >> headers.txt
        fi
        if ! bcftools view -h annotated.vcf.gz | grep -q '##INFO=<ID=SVTYPE'; then
            echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">' >> headers.txt
        fi

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' annotated.vcf.gz | \
        awk '{
            ref_len = length($3)
            alt_len = length($4)
            if (ref_len == alt_len) {
                svtype = "SNV"
            } else if (ref_len > alt_len) {
                svtype = "DEL"
            } else {
                svtype = "INS"
            }
            print $1"\t"$2"\t"$3"\t"$4"\t"svtype
        }' | bgzip -c > svtype_annot.txt.gz
        tabix -s1 -b2 -e2 svtype_annot.txt.gz

        header_flag=""
        if [ -s headers.txt ]; then
            header_flag="-h headers.txt"
        fi

        bcftools annotate -a svtype_annot.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/SVTYPE \
            $header_flag \
            annotated.vcf.gz -Oz -o ~{prefix}.vcf.gz

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

