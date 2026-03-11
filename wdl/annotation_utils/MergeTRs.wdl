version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tr_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetBase {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_base
        }

        call Helpers.SubsetVcfToSamples as SubsetSamplesBase {
            input:
                vcf = SubsetContigBase.subset_vcf,
                vcf_idx = SubsetContigBase.subset_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.subset_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_base
        }

        scatter (i in range(length(tr_vcfs))) {
            call Helpers.SubsetVcfToContig as SubsetContigTr {
                input:
                    vcf = tr_vcfs[i],
                    vcf_idx = tr_vcf_idxs[i],
                    contig = contig,
                    extra_args = "--min-ac 1",
                    prefix = "~{prefix}.~{contig}.tr~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig_tr
            }

            call Helpers.SubsetVcfToSamples as SubsetSamplesTr {
                input:
                    vcf = tr_vcfs[i],
                    vcf_idx = tr_vcf_idxs[i],
                    samples = sample_ids,
                    prefix = "~{prefix}.~{contig}.tr~{i}.subset_samples",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_samples_tr
            }

            call Helpers.CheckSampleConsistency {
                input:
                    vcfs = [SubsetSamplesBase.subset_vcf, SubsetSamplesTr.subset_vcf],
                    vcf_idxs = [SubsetSamplesBase.subset_vcf_idx, SubsetSamplesTr.subset_vcf_idx],
                    sample_ids = sample_ids,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_check_sample_consistency
            }

            call Helpers.SetMissingFiltersToPass as SetTrFiltersToPass {
                input:
                    vcf = SubsetSamplesTr.subset_vcf,
                    vcf_idx = SubsetSamplesTr.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.refiltered",
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

            call DeduplicateOverlappingVariants {
                input:
                    vcf = TagTRVcf.tagged_vcf,
                    vcf_idx = TagTRVcf.tagged_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.dedup",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_deduplicate_trs
            }
        }

        call PriorityMergeTRVcfs {
            input:
                vcfs = DeduplicateOverlappingVariants.dedup_vcf,
                vcf_idxs = DeduplicateOverlappingVariants.dedup_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_priority_merge
        }

        call SetTRVariantIds {
            input:
                vcf = PriorityMergeTRVcfs.merged_vcf,
                vcf_idx = PriorityMergeTRVcfs.merged_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr.ids",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_tr_ids
        }

        call AnnotateVcfWithTRs {
            input:
                vcf = SubsetSamplesBase.subset_vcf,
                vcf_idx = SubsetSamplesBase.subset_vcf_idx,
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
            allow_overlaps = false,
            naive = true,
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

task DeduplicateOverlappingVariants {
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

task PriorityMergeTRVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCF_FILES=(~{sep=' ' vcfs})

        cp "${VCF_FILES[0]}" merged.vcf.gz

        tabix -p vcf merged.vcf.gz

        for vcf_file in "${VCF_FILES[@]:1}"; do
            # BED of variants currently in the merged set
            bcftools query \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
                merged.vcf.gz \
                | awk 'BEGIN{OFS="\t"} {print $1, $2, $2+length($4)}' > merged.bed

            # BED of candidate variants
            bcftools query \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
                "$vcf_file" \
                | awk 'BEGIN{OFS="\t"} {print $1, $2, $2+length($4), $3}' > candidates.bed

            # Keep only candidates with no overlap in the merged set
            bedtools intersect \
                -v \
                -a candidates.bed \
                -b merged.bed \
                | cut -f4 > ids_to_add.txt

            if [ -s ids_to_add.txt ]; then
                bcftools view \
                    -i "ID=@ids_to_add.txt" \
                    -Oz -o to_add.vcf.gz \
                    "$vcf_file"
                
                tabix -p vcf to_add.vcf.gz

                bcftools concat \
                    --allow-overlaps \
                    merged.vcf.gz to_add.vcf.gz \
                    | bcftools sort -Oz -o new_merged.vcf.gz
                
                tabix -p vcf new_merged.vcf.gz

                mv new_merged.vcf.gz merged.vcf.gz

                mv new_merged.vcf.gz.tbi merged.vcf.gz.tbi

                rm -f to_add.vcf.gz to_add.vcf.gz.tbi
            fi

            rm -f merged.bed candidates.bed ids_to_add.txt
        done

        mv merged.vcf.gz ~{prefix}.vcf.gz
        mv merged.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcfs, "GB")) + 10,
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
