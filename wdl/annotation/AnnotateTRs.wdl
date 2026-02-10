version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        String tr_caller

        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tr_vcf
        RuntimeAttr? runtime_attr_set_missing_filters
        RuntimeAttr? runtime_attr_set_tr_ids
        RuntimeAttr? runtime_attr_annotate_trs
        RuntimeAttr? runtime_attr_concat_vcf
    }

    call CheckSamplesMatch {
        input:
            vcf = vcf,
            tr_vcf = tr_vcf,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
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

        call Helpers.SubsetVcfToContig as SubsetTrVcf {
            input:
                vcf = tr_vcf,
                vcf_idx = tr_vcf_idx,
                contig = contig,
                extra_args = "--min-ac 1",
                prefix = "~{prefix}.~{contig}.tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tr_vcf
        }

        call Helpers.SetMissingFiltersToPass {
            input:
                vcf = SubsetTrVcf.subset_vcf,
                vcf_idx = SubsetTrVcf.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr.pass",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_missing_filters
        }

        call SetTRVariantIds {
            input:
                vcf = SetMissingFiltersToPass.filtered_vcf,
                vcf_idx = SetMissingFiltersToPass.filtered_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr.ids",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_tr_ids
        }

        call AnnotateTRVariants {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                tr_vcf = SetTRVariantIds.renamed_vcf,
                tr_vcf_idx = SetTRVariantIds.renamed_vcf_idx,
                tr_caller = tr_caller,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                check_passed = CheckSamplesMatch.check_passed,
                runtime_attr_override = runtime_attr_annotate_trs
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateTRVariants.annotated_vcf,
            vcf_idxs = AnnotateTRVariants.annotated_vcf_idx,
            prefix = prefix + ".annotated_trs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File tr_annotated_vcf = ConcatVcfs.concat_vcf
        File tr_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task CheckSamplesMatch {
    input {
        File vcf
        File tr_vcf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} | sort > samples_vcf.txt
        bcftools query -l ~{tr_vcf} | sort > samples_tr_vcf.txt

        if ! diff samples_vcf.txt samples_tr_vcf.txt > /dev/null; then
            echo "ERROR: Sample lists do not match between VCF and TR VCF" >&2
            echo "VCF samples:" >&2
            cat samples_vcf.txt >&2
            echo "TR VCF samples:" >&2
            cat samples_tr_vcf.txt >&2
            exit 1
        fi

        echo "Sample lists match"
        echo "true" > check_passed.txt
    >>>

    output {
        Boolean check_passed = read_boolean("check_passed.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 5,
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

# First pass: count ID occurrences
vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
    id_counts[new_id] += 1
vcf_in.close()

# Second pass: assign IDs with suffixes if needed
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

task AnnotateTRVariants {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        String tr_caller
        String prefix
        String docker
        Boolean check_passed
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # TR VCF processing
        bcftools query \
            -l \
            ~{vcf} \
            > samples.txt
        
        bcftools view \
            -S samples.txt \
            -Oz -o tr_reordered.vcf.gz \
            ~{tr_vcf}
        
        tabix -p vcf tr_reordered.vcf.gz

        touch new_headers.txt
        if ! bcftools view -h tr_reordered.vcf.gz | grep -q '##INFO=<ID=allele_type'; then
            echo '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">' >> new_headers.txt
        fi
        if ! bcftools view -h tr_reordered.vcf.gz | grep -q '##INFO=<ID=SOURCE'; then
            echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of variant call">' >> new_headers.txt
        fi
        
        bcftools annotate \
            -h new_headers.txt \
            -Oz -o tr_header_added.vcf.gz \
            tr_reordered.vcf.gz

        bcftools view tr_header_added.vcf.gz \
            | awk -v source="~{tr_caller}" 'BEGIN{OFS="\t"} /^#/ {print; next} {
                $8 = $8 ";allele_type=trv;SOURCE=" source
                print
            }' | bgzip -c > tr_tagged.vcf.gz
        
        tabix -p vcf tr_tagged.vcf.gz

        # Extract overlaps
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            tr_tagged.vcf.gz \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                print $1, $2, end, $3, $4, $5
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

        # Annotate overlaps
        awk -v filter="~{tr_caller}_OVERLAPPED" 'BEGIN{OFS="\t"} {
            print $1, $2, $5, $6, $4, filter, $10
        }' overlaps.bed | sort -k1,1 -k2,2n | bgzip -c > annotations.tsv.gz

        tabix -s 1 -b 2 -e 2 annotations.tsv.gz

        bcftools view \
            -h ~{vcf} \
        | grep -v 'ID=AL,' > vcf_no_al_header.txt

        bcftools view \
            -h ~{tr_vcf} \
        | grep 'ID=AL,' > al_header.txt || true

        bcftools reheader \
            -h vcf_no_al_header.txt \
            ~{vcf} \
            > vcf_stripped.vcf.gz

        cat <<EOF > new_header.txt
##INFO=<ID=TRID,Number=1,Type=String,Description="ID of enveloping tandem repeat">
##FILTER=<ID=~{tr_caller}_OVERLAPPED,Description="Variant enveloped by ~{tr_caller}">
EOF
        
        cat new_header.txt al_header.txt > merged_headers.txt

        bcftools annotate \
            -a annotations.tsv.gz \
            -c CHROM,POS,REF,ALT,~ID,=FILTER,INFO/TRID \
            -h merged_headers.txt \
            -Oz -o vcf_annotated.vcf.gz \
            vcf_stripped.vcf.gz

        tabix -p vcf vcf_annotated.vcf.gz

        bcftools concat \
            --allow-overlaps \
            vcf_annotated.vcf.gz \
            tr_tagged.vcf.gz \
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
        disk_gb: 5 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 20,
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
