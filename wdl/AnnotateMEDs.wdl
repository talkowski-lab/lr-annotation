version 1.0

import "general/Structs.wdl"
import "general/Helpers.wdl"

workflow AnnotateMEDs {
    input {
        File vcf
        File vcf_idx
        File contig_list
        File med_catalog

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

    Array[String] contigs = read_lines(contig_list)

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetContig {
            input:
                vcf = vcf,
                vcf_index = vcf_idx,
                contig = contig,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call ExtractDeletionsToBed as ExtractDEL {
            input:
                vcf = SubsetContig.subset_vcf,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bedtools
        }

        call BedtoolsIntersect as IntersectMED {
            input:
                bed_a = ExtractDEL.del_bed,
                bed_b = med_catalog,
                prefix = "~{prefix}.~{contig}",
                reciprocal_overlap = reciprocal_overlap,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bedtools
        }

        call AnnotateVcfWithMED as AnnotateContig {
            input:
                vcf = SubsetContig.subset_vcf,
                intersect_bed = IntersectMED.intersect_bed,
                med_catalog = med_catalog,
                prefix = "~{prefix}.~{contig}",
                size_similarity = size_similarity,
                reciprocal_overlap = reciprocal_overlap,
                breakpoint_window = breakpoint_window,
                sequence_similarity = sequence_similarity,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs as ConcatVCF {
        input:
            vcfs = AnnotateContig.annotated_vcf,
            vcfs_idx = AnnotateContig.annotated_vcf_idx,
            outfile_prefix = prefix + ".med_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File med_vcf = ConcatVCF.concat_vcf
        File med_vcf_idx = ConcatVCF.concat_vcf_idx
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
            bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' > ~{prefix}.del.bed
    >>>

    output {
        File del_bed = "~{prefix}.del.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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

task BedtoolsIntersect {
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
        mem_gb: 2,
        disk_gb: 10,
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

task AnnotateVcfWithMED {
    input {
        File vcf
        File intersect_bed
        File med_catalog
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

        if bcftools view -h ~{vcf} | grep -q "ID=ME_TYPE"; then
            touch me_header.txt
        else
            echo '##INFO=<ID=ME_TYPE,Number=1,Type=String,Description="Type of mobile element">' > me_header.txt
        fi

        python3 <<'EOF'
import sys
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
    # A (VCF): 0=chrom, 1=start, 2=end, 3=ID
    # B (MED): 4=chrom, 5=start, 6=end, 7=ID, 8=seq, 9=designation
    del_chrom = row[0]
    del_start = int(row[1])
    del_end = int(row[2])
    del_id = row[3]

    med_start = int(row[5])
    med_end = int(row[6])
    designation = row[9]
    
    me_type = get_me_type(designation)
    
    if me_type and passes_criteria(del_start, del_end, med_start, med_end, size_similarity, reciprocal_overlap, breakpoint_window):
        annotations.append([del_chrom, del_start + 1, del_id, me_type])

out_df = pd.DataFrame(annotations, columns=['CHROM', 'POS', 'ID', 'ME_TYPE'])
out_df.drop_duplicates(subset=['ID'], inplace=True)
out_df.to_csv("annotations.tsv", sep='\t', index=False, header=False)

EOF

        sort -k1,1 -k2,2n annotations.tsv | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz

        bcftools annotate \
            -a annotations.tsv.gz \
            -h me_header.txt \
            -c CHROM,POS,ID,INFO/ME_TYPE \
            -O z \
            -o ~{prefix}.med_annotated.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.med_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.med_annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.med_annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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