version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateRegion {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masker_bed
        Float min_coverage

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_label_regions
        RuntimeAttr? runtime_attr_set_unique
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call LabelVariantRegions as LabelSimpleRepeats {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                region_bed = simple_repeats_bed,
                region_name = "SR",
                min_coverage = min_coverage,
                prefix = "~{prefix}.~{contig}.simple_repeats",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_label_regions
        }

        call LabelVariantRegions as LabelSegDups {
            input:
                vcf = LabelSimpleRepeats.labeled_vcf,
                vcf_idx = LabelSimpleRepeats.labeled_vcf_idx,
                region_bed = seg_dup_bed,
                region_name = "SD",
                min_coverage = min_coverage,
                prefix = "~{prefix}.~{contig}.seg_dups",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_label_regions
        }

        call LabelVariantRegions as LabelRepeatMasker {
            input:
                vcf = LabelSegDups.labeled_vcf,
                vcf_idx = LabelSegDups.labeled_vcf_idx,
                region_bed = repeat_masker_bed,
                region_name = "RM",
                min_coverage = min_coverage,
                prefix = "~{prefix}.~{contig}.repeat_masker",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_label_regions
        }

        call SetUniqueRegion {
            input:
                vcf = LabelRepeatMasker.labeled_vcf,
                vcf_idx = LabelRepeatMasker.labeled_vcf_idx,
                prefix = "~{prefix}.~{contig}.region_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_unique
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = SetUniqueRegion.annotated_vcf,
            vcf_idxs = SetUniqueRegion.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.region_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File region_annotated_vcf = ConcatVcfs.concat_vcf
        File region_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task LabelVariantRegions {
    input {
        File vcf
        File vcf_idx
        File region_bed
        String region_name
        Float min_coverage
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Add REGION header if not already present
        touch header.lines
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=REGION'; then
            echo '##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region context (e.g. SIMPLE_REPEAT, SEG_DUP, REPEAT_MASKER). UNIQUE if no region overlaps by the minimum coverage threshold.">' \
                > header.lines
        fi
        bcftools annotate \
            -h header.lines \
            -Oz -o with_header.vcf.gz \
            ~{vcf}
        tabix -p vcf with_header.vcf.gz

        # Extract variants
        bcftools query \
            -f '%CHROM\t%POS0\t%END\t%ID\t%REF\t%ALT\n' \
            -e 'INFO/REGION!="."' \
            with_header.vcf.gz \
        | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"($3-$2)"\t"$5"\t"$6}' \
        | sort -k1,1 -k2,2n > variants.bed

        # Intersect variants with regions
        bedtools intersect \
            -a variants.bed \
            -b ~{region_bed} \
            -wo \
        | awk -v min_cov=~{min_coverage} -v region="~{region_name}" \
            '{ vid=$4; vlen=($5>0?$5:1); overlap=$NF
               if (overlap/vlen >= min_cov && !seen[vid]++)
                   print $1"\t"($2+1)"\t"$6"\t"$7"\t"vid"\t"region }' \
        | sort -k1,1 -k2,2n \
        | bgzip -c > annot.txt.gz

        tabix -s1 -b2 -e2 annot.txt.gz

        # Annotate flagged variants
        bcftools annotate \
            -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/REGION \
            with_header.vcf.gz \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File labeled_vcf = "~{prefix}.vcf.gz"
        File labeled_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(region_bed, "GB")) + 10,
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

task SetUniqueRegion {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Annotate any variants without INFO/REGION as UNIQUE
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            -e 'INFO/REGION!="."' \
            ~{vcf} \
        | awk '{print $0"\tUS"}' \
        | bgzip -c > unique_annot.txt.gz

        tabix -s1 -b2 -e2 unique_annot.txt.gz

        bcftools annotate \
            -a unique_annot.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/REGION \
            ~{vcf} \
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
