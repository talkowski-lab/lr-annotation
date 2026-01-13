version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateMEDs {
    input {
        File vcf
        File vcf_idx
        File med_catalog
        Array[String] contigs

        String prefix
        String utils_docker

        Float size_similarity = 0.5
        Float reciprocal_overlap = 0.5
        Int breakpoint_window = 5000
        Float sequence_similarity = 0.5

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_bedtools
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_index = vcf_idx,
                contig = contig,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call ExtractDeletionsToBed {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bedtools
        }

        call IntersectMED {
            input:
                bed_a = ExtractDeletionsToBed.del_bed,
                bed_b = med_catalog,
                prefix = "~{prefix}.~{contig}",
                reciprocal_overlap = reciprocal_overlap,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bedtools
        }

        call GenerateMedAnnotationTable {
            input:
                intersect_bed = IntersectMED.intersect_bed,
                prefix = "~{prefix}.~{contig}",
                size_similarity = size_similarity,
                reciprocal_overlap = reciprocal_overlap,
                breakpoint_window = breakpoint_window,
                sequence_similarity = sequence_similarity,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatTsvs as MergeAnnotations {
        input:
            tsvs = GenerateMedAnnotationTable.annotations_tsv,
            prefix = prefix + ".med_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File annotations_tsv_meds = MergeAnnotations.concatenated_tsv
    }
}

task ExtractDeletionsToBed {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -i 'SVTYPE=="DEL"' ~{vcf} | \
            bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%ID\n' > ~{prefix}.del.bed
    >>>

    output {
        File del_bed = "~{prefix}.del.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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

task IntersectMED {
    input {
        File bed_a
        File bed_b
        String prefix
        Float reciprocal_overlap
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bedtools intersect \
            -wa -wb \
            -r -f ~{reciprocal_overlap} \
            -a ~{bed_a} \
            -b ~{bed_b} \
            > ~{prefix}.intersect.bed
    >>>

    output {
        File intersect_bed = "~{prefix}.intersect.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_a, "GB") + size(bed_b, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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

task GenerateMedAnnotationTable {
    input {
        File intersect_bed
        String prefix
        Float size_similarity
        Float reciprocal_overlap
        Int breakpoint_window
        Float sequence_similarity
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import pandas as pd

def get_me_type(designation):
    d = str(designation)
    if d.startswith('SINE') or d.startswith('Alu'):
        return 'ALU'
    elif d.startswith('LINE') or d.startswith('L1'):
        return 'LINE'
    elif d.startswith('Retroposon') or d.startswith('SVA'):
        return 'SVA'
    return None

def passes_criteria(del_start, del_end, med_start, med_end, size_similarity, reciprocal_overlap, breakpoint_window):
    del_len = abs(del_end - del_start)
    med_len = abs(med_end - med_start)
    
    if del_len == 0 or med_len == 0:
        return False
        
    size_sim = min(del_len, med_len) / max(del_len, med_len)
    overlap = max(0, min(del_end, med_end) - max(del_start, med_start))
    rec_ovl = overlap / min(del_len, med_len)
    
    if size_sim < size_similarity:
        return False
    if rec_ovl < reciprocal_overlap:
        return False
    if abs(del_start - med_start) > breakpoint_window and abs(del_end - med_end) > breakpoint_window:
        return False

    return True

bed = pd.read_csv("~{intersect_bed}", sep='\t', header=None)

size_similarity = float(~{size_similarity})
reciprocal_overlap = float(~{reciprocal_overlap})
breakpoint_window = int(~{breakpoint_window})

annotations = []

for _, row in bed.iterrows():
    # bedtools intersect -wa -wb outputs:
    # A (del.bed): 0=chrom, 1=pos, 2=end, 3=ref, 4=alt, 5=id
    # B (MED): 6=chrom, 7=start, 8=end, 9=ID, 10=seq, 11=designation
    del_chrom = row[0]
    del_start = int(row[1])
    del_end = int(row[2])
    del_ref = row[3]
    del_alt = row[4]
    del_id = row[5]

    med_start = int(row[7])
    med_end = int(row[8])
    designation = row[11]
    
    me_type = get_me_type(designation)
    
    if not me_type:
        continue

    if passes_criteria(del_start, del_end, med_start, med_end, size_similarity, reciprocal_overlap, breakpoint_window):
        annotations.append([
            del_chrom,
            del_start,
            del_ref,
            del_alt,
            del_id,
            me_type
        ])

out_df = pd.DataFrame(annotations, columns=['CHROM', 'POS', 'REF', 'ALT', 'ID', 'ME_TYPE'])
out_df.drop_duplicates(subset=['ID'], inplace=True)
out_df.sort_values(['CHROM', 'POS'], inplace=True)
out_df.to_csv("~{prefix}.annotations.tsv", sep='\t', index=False, header=False)

EOF
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(intersect_bed, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
