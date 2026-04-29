version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateRegion {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masker_bed
        Float min_coverage

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_label_regions
        RuntimeAttr? runtime_attr_set_unique
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcf
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords as ShardVcf {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.region_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcf.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcf.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call LabelVariantRegions as LabelSimpleRepeats {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    region_bed = simple_repeats_bed,
                    region_name = "SR",
                    min_coverage = min_coverage,
                    prefix = "~{prefix}.~{contig}.simple_repeats.shard_~{i}",
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
                    prefix = "~{prefix}.~{contig}.seg_dups.shard_~{i}",
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
                    prefix = "~{prefix}.~{contig}.repeat_masker.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_label_regions
            }

            call SetUniqueRegion {
                input:
                    vcf = LabelRepeatMasker.labeled_vcf,
                    vcf_idx = LabelRepeatMasker.labeled_vcf_idx,
                    prefix = "~{prefix}.~{contig}.region_annotated.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_set_unique
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = SetUniqueRegion.annotations_tsv,
                    sort_output = true,
                    prefix = "~{prefix}.~{contig}.region_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, SetUniqueRegion.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.region_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_region = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
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
            echo '##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region context.">' \
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

        # Extract all variants with their REGION, defaulting to US for unset
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/REGION\n' \
            ~{vcf} \
        | awk 'BEGIN{OFS="\t"} {if ($6 == ".") $6 = "US"; print}' \
            > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
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
