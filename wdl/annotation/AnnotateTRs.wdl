version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateTRs {
    input {
        File vcf
        File vcf_idx
        Array[File] tr_vcfs
        Array[File] tr_vcf_idxs
        Array[String] contigs
        String prefix

        Array[String] tr_callers

        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tr_vcf
        RuntimeAttr? runtime_attr_set_missing_filters
        RuntimeAttr? runtime_attr_tag_tr_vcf
        RuntimeAttr? runtime_attr_deduplicate_trs
        RuntimeAttr? runtime_attr_concat_tr_vcfs
        RuntimeAttr? runtime_attr_set_tr_ids
        RuntimeAttr? runtime_attr_annotate_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (i in range(length(tr_vcfs))) {
        call Helpers.CheckSampleMatch {
            input:
                vcf_a = vcf,
                vcf_a_idx = vcf_idx,
                vcf_b = tr_vcfs[i],
                vcf_b_idx = tr_vcf_idxs[i],
                docker = utils_docker,
                runtime_attr_override = runtime_attr_check_samples
        }
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        scatter (i in range(length(tr_vcfs))) {
            call Helpers.SubsetVcfToContig as SubsetTrVcf {
                input:
                    vcf = tr_vcfs[i],
                    vcf_idx = tr_vcf_idxs[i],
                    contig = contig,
                    extra_args = "--min-ac 1",
                    prefix = "~{prefix}.~{contig}.tr~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tr_vcf
            }

            call Helpers.SetMissingFiltersToPass as SetTrFiltersToPass {
                input:
                    vcf = SubsetTrVcf.subset_vcf,
                    vcf_idx = SubsetTrVcf.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.pass",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_set_missing_filters
            }

            call TagTRVcf {
                input:
                    vcf = SetTrFiltersToPass.filtered_vcf,
                    vcf_idx = SetTrFiltersToPass.filtered_vcf_idx,
                    tr_caller = tr_callers[i],
                    prefix = "~{prefix}.~{contig}.tr~{i}.tagged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_tag_tr_vcf
            }
        }

        call Helpers.ConcatVcfs as ConcatTrVcfs {
            input:
                vcfs = TagTRVcf.tagged_vcf,
                vcf_idxs = TagTRVcf.tagged_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_concat",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_tr_vcfs
        }

        call DeduplicateEnvelopedVariants {
            input:
                vcf = ConcatTrVcfs.concat_vcf,
                vcf_idx = ConcatTrVcfs.concat_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_dedup",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_deduplicate_trs
        }

        call SetTRVariantIds {
            input:
                vcf = DeduplicateEnvelopedVariants.dedup_vcf,
                vcf_idx = DeduplicateEnvelopedVariants.dedup_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr.ids",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_tr_ids
        }

        call AnnotateVcfWithTRs {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                tr_vcf = SetTRVariantIds.renamed_vcf,
                tr_vcf_idx = SetTRVariantIds.renamed_vcf_idx,
                tr_callers = tr_callers,
                prefix = "~{prefix}.~{contig}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateVcfWithTRs.annotated_vcf,
            vcf_idxs = AnnotateVcfWithTRs.annotated_vcf_idx,
            prefix = "~{prefix}.annotated_trs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File tr_annotated_vcf = ConcatVcfs.concat_vcf
        File tr_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task TagTRVcf {
    input {
        File vcf
        File vcf_idx
        String tr_caller
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -l \
            ~{vcf} \
            > samples.txt
        
        bcftools view \
            -S samples.txt \
            -Oz -o tr_reordered.vcf.gz \
            ~{vcf}
        
        tabix -p vcf tr_reordered.vcf.gz

        touch new_headers.txt
        if ! bcftools view -h tr_reordered.vcf.gz | grep -q '##INFO=<ID=allele_type'; then
            echo '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">' >> new_headers.txt
        fi
        if ! bcftools view -h tr_reordered.vcf.gz | grep -q '##INFO=<ID=SOURCE'; then
            echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of variant call">' >> new_headers.txt
        fi
        if ! bcftools view -h tr_reordered.vcf.gz | grep -q '##INFO=<ID=TR_OVERLAPPED'; then
            echo '##INFO=<ID=TR_OVERLAPPED,Number=0,Type=Flag,Description="Variant enveloped by tandem repeat">' >> new_headers.txt
        fi
        
        bcftools annotate \
            -h new_headers.txt \
            -Oz -o tr_header_added.vcf.gz \
            tr_reordered.vcf.gz

        rm -f tr_reordered.vcf.gz new_headers.txt

        bcftools view tr_header_added.vcf.gz \
            | awk -v source="~{tr_caller}" 'BEGIN{OFS="\t"} /^#/ {print; next} {
                $3 = $1 "-" $2 "-" source "-" length($4)
                $8 = $8 ";allele_type=trv;SOURCE=" source
                print
            }' \
            | bgzip -c > ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
        
        rm -f tr_header_added.vcf.gz
    >>>

    output {
        File tagged_vcf = "~{prefix}.vcf.gz"
        File tagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 5,
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

task DeduplicateEnvelopedVariants {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                ref_len = length($4)
                print $1, $2, end, $3, ref_len
            }' > variants.bed

        bedtools intersect \
            -f 1.0 \
            -wa \
            -wb \
            -a variants.bed \
            -b variants.bed \
            | awk 'BEGIN{OFS="\t"} $4 != $9 {
                if ($5 < $10) {
                    print $4
                } else if ($5 == $10 && $4 > $9) {
                    print $4
                }
            }' \
            | sort -u > ids_to_remove.txt

        if [ -s ids_to_remove.txt ]; then
            bcftools view \
                -e "ID=@ids_to_remove.txt" \
                -Oz -o ~{prefix}.vcf.gz \
                ~{vcf}
        else
            cp ~{vcf} ~{prefix}.vcf.gz
        fi
        
        tabix -p vcf ~{prefix}.vcf.gz
        
        rm -f variants.bed ids_to_remove.txt
    >>>

    output {
        File dedup_vcf = "~{prefix}.vcf.gz"
        File dedup_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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

task SetTRVariantIds {
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
from collections import defaultdict

vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
    id_counts[new_id] += 1
vcf_in.close()

vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
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
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task AnnotateVcfWithTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] tr_callers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SOURCE\n' \
            ~{tr_vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                print $1, $2, end, $3, $4, $5, $6
            }' > tr.bed
        
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                print $1, $2, end, $3, $4, $5
            }' > vcf.bed

        bedtools intersect \
            -f 1.0 \
            -wa \
            -wb \
            -a vcf.bed \
            -b tr.bed \
            > overlaps.bed

        rm -f tr.bed vcf.bed

        awk 'BEGIN{OFS="\t"} {
            print $1, $2, $5, $6, $4, "1", $10
        }' overlaps.bed \
            | sort -k1,1 -k2,2n \
            | bgzip -c > annotations.tsv.gz

        tabix -s 1 -b 2 -e 2 annotations.tsv.gz
        
        rm -f overlaps.bed

        bcftools view \
            -h \
            ~{vcf} \
        | grep -v 'ID=AL,' > vcf_no_al_header.txt

        bcftools view \
            -h \
            ~{tr_vcf} \
        | grep 'ID=AL,' > al_header.txt || true

        bcftools reheader \
            -h vcf_no_al_header.txt \
            ~{vcf} \
            > vcf_stripped.vcf.gz

        cat <<EOF > new_header.txt
##INFO=<ID=TR_OVERLAPPED,Number=0,Type=Flag,Description="Variant enveloped by tandem repeat">
##INFO=<ID=TRID,Number=1,Type=String,Description="ID of enveloping tandem repeat">
EOF
        
        cat new_header.txt al_header.txt > merged_headers.txt

        bcftools annotate \
            -a annotations.tsv.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/TR_OVERLAPPED,INFO/TRID \
            -h merged_headers.txt \
            -Oz -o vcf_annotated.vcf.gz \
            vcf_stripped.vcf.gz

        tabix -p vcf vcf_annotated.vcf.gz
        
        rm -f annotations.tsv.gz annotations.tsv.gz.tbi vcf_no_al_header.txt al_header.txt new_header.txt merged_headers.txt vcf_stripped.vcf.gz

        bcftools concat \
            --allow-overlaps \
            -Oz -o ~{prefix}.vcf.gz \
            vcf_annotated.vcf.gz \
            ~{tr_vcf}

        tabix -p vcf ~{prefix}.vcf.gz
        
        rm -f vcf_annotated.vcf.gz vcf_annotated.vcf.gz.tbi
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 20,
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
