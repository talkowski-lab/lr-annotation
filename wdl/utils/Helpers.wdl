version 1.0

import "Structs.wdl"

task AddFilter {
    input {
        File vcf
        File vcf_idx
        String filter_name
        String filter_description
        String filter_expression
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -h \
            ~{vcf} \
        | grep "^##" > header.txt
        
        echo '##FILTER=<ID=~{filter_name},Description="~{filter_description}">' >> header.txt
        
        bcftools view \
            -h \
            ~{vcf} \
        | grep "^#CHROM" >> header.txt
        
        bcftools reheader \
            -h header.txt \
            ~{vcf} \
        | bcftools filter \
            --mode + \
            -s ~{filter_name} \
            -e '~{filter_expression}' \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File flagged_vcf = "~{prefix}.vcf.gz"
        File flagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task AddInfo {
    input {
        File vcf
        File vcf_idx
        String tag_id
        String tag_value
        String tag_description
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=~{tag_id},Number=1,Type=String,Description="~{tag_description}">' > header.lines

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t~{tag_value}\n' ~{vcf} | \
        bgzip -c > annotations.txt.gz
        
        tabix -s1 -b2 -e2 annotations.txt.gz

        bcftools annotate -h header.lines -a annotations.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/~{tag_id} \
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

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String svtype
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) | sed -e "s/#//g" > ~{svtype}.bed
        bedtools closest -wo -a <(sort -k1,1 -k2,2n ~{bed_a}) -b <(sort -k1,1 -k2,2n ~{bed_b}) >> ~{svtype}.bed
    >>>

    output {
        File output_bed = "~{svtype}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_a, "GB") + size(bed_b, "GB")) + 5,
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

task ConcatAlignedTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
import csv

input_files = "~{sep=',' tsvs}".split(',')
output_filename = "aligned_unsorted.tsv"
header_filename = "~{prefix}.header.txt"
fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]

all_keys = set()
for f in input_files:
    with open(f, 'r') as fh:
        line = fh.readline().strip()
        if not line: continue
        parts = line.split('\t')
        if len(parts) > 5:
            keys = parts[5:]
            all_keys.update(keys)

sorted_keys = sorted(list(all_keys))
master_header = fixed_cols + sorted_keys

with open(header_filename, 'w') as hout:
    for k in sorted_keys:
        hout.write(k + "\n")

with open(output_filename, 'w') as out:
    out.write("\t".join(master_header) + "\n")
    
    for f in input_files:
        with open(f, 'r') as fh:
            header_line = fh.readline().strip()
            if not header_line: 
                continue
            
            file_cols = header_line.split('\t')
            col_map = {name: i for i, name in enumerate(file_cols)}
            
            for line in fh:
                parts = line.strip().split('\t')
                if not parts: continue
                
                out_row = []
                for target_col in master_header:
                    if target_col in col_map:
                        try:
                            val = parts[col_map[target_col]]
                            out_row.append(val)
                        except IndexError:
                            out_row.append(".")
                    else:
                        out_row.append(".")
                
                out.write("\t".join(out_row) + "\n")
CODE
    
        tail -n +2 aligned_unsorted.tsv | sort -k1,1 -k2,2n > ~{prefix}.tsv
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
        File merged_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
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

task ConcatTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        Boolean preserve_header = false
        Boolean skip_sort = false
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ "~{preserve_header}" == "true" ]; then
            head -n 1 ~{tsvs[0]} > ~{prefix}.tsv
            tail -n +2 -q ~{sep=' ' tsvs} >> ~{prefix}.tsv
        elif [ "~{skip_sort}" == "true" ]; then
            cat ~{sep=' ' tsvs} > ~{prefix}.tsv
        else
            cat ~{sep=' ' tsvs} | sort -k1,1 -k2,2n > ~{prefix}.tsv
        fi
    >>>

    output {
        File concatenated_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
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

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcfs_idx
        Boolean merge_sort = true
        String prefix = "concat"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String merge_flag = if merge_sort then "--allow-overlaps" else ""

    command <<<
            set -euo pipefail
            
            VCFS_FILE="~{write_lines(vcfs)}"

            bcftools concat \
                ~{merge_flag} \
                --file-list ${VCFS_FILE} \
                -Oz -o "~{prefix}.vcf.gz"
            
            tabix -p vcf -f "~{prefix}.vcf.gz"
    >>>

    output {
            File concat_vcf = "~{prefix}.vcf.gz"
            File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
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

task ConvertToSymbolic {
    input {
        File vcf
        File vcf_idx
        String prefix
        Boolean drop_genotypes = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -f '%INFO/SVTYPE\n' ~{vcf} | sort -u > svtypes.txt

        if [ "~{drop_genotypes}" == "true" ]; then
            bcftools view \
                -G \
                ~{vcf} \
            | python3 /opt/gnomad-lr/scripts/helpers/symbalts.py \
                - \
                svtypes.txt \
            | bcftools view \
                -Oz -o ~{prefix}.vcf.gz
        else
            python3 /opt/gnomad-lr/scripts/helpers/symbalts.py \
                ~{vcf} \
                svtypes.txt \
            | bcftools view \
                -Oz -o ~{prefix}.vcf.gz
        fi
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task ExtractVcfAnnotations {
    input {
        File vcf
        File vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        Boolean add_header_row = false
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import sys

vcf = pysam.VariantFile("~{vcf}")
orig = pysam.VariantFile("~{original_vcf}")

new_keys = sorted(list(set(vcf.header.info.keys()) - set(orig.header.info.keys())))

with open("~{prefix}.header.txt", "w") as out:
    for k in new_keys:
        out.write(k + "\n")

with open("~{prefix}.annotations.tsv", "w") as out:
    
    if "~{add_header_row}" == "true":
        fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]
        full_header = fixed_cols + new_keys
        out.write("\t".join(full_header) + "\n")

    for record in vcf:
        alts = ",".join(record.alts) if record.alts else "."
        rid = record.id if record.id else "."
        row = [record.chrom, str(record.pos), record.ref, alts, rid]

        for k in new_keys:
            if k in record.info:
                val = record.info[k]
                if isinstance(val, bool):
                    row.append("1" if val else "0")
                elif isinstance(val, (list, tuple)):
                    row.append(",".join(map(str, val)))
                else:
                    row.append(str(val))
            else:
                row.append(".")
        
        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
        File annotations_header = "~{prefix}.header.txt"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(original_vcf, "GB")) + 5,
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

task FinalizeToDir {
    input {
        Array[File] files
        String outdir
        File? keyfile
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    command <<<
        set -euo pipefail

        cat ~{write_lines(files)} | gsutil -m cp -I "~{gcs_output_dir}"
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(files, "GB")) + 5,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FinalizeToFile {
    input {
        File file
        String outdir
        String? name
        File? keyfile
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(file, "GB")) + 5,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GetHailMTSize {
    input {
        String mt_uri
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        tot_size=$(gsutil -m du -sh ~{mt_uri} | awk -F '    ' '{ print $1 }')

        python3 <<CODE > mt_size.txt
import sys

size = "$tot_size".split()[0]
unit = "$tot_size".split()[1]

def convert_to_gib(size, unit):
    size_dict = {"KiB": 2**10, "MiB": 2**20, "GB": 2**30, "TiB": 2**40}
    return float(size) * size_dict[unit] / size_dict["GB"]

size_in_gib = convert_to_gib(size, unit)
print(size_in_gib)
CODE
    >>>

    output {
        Float mt_size = read_lines('mt_size.txt')[0]
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: 25,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
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

task MergeHeaderLines {
    input {
        Array[File] header_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=' ' header_files} \
            | sort -u > ~{prefix}.txt
    >>>

    output {
        File merged_header = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task RenameVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String id_format = "%CHROM-%POS-%REF-%ALT"
        Boolean strip_chr = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools annotate --set-id '~{id_format}' ~{vcf} -Oz -o temp_renamed.vcf.gz
        
        if [ "~{strip_chr}" == "true" ]; then
            bcftools view -h temp_renamed.vcf.gz > header.txt
            bcftools view -H temp_renamed.vcf.gz \
                | awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $3); print}' \
                | cat header.txt - \
                | bgzip -c > ~{prefix}.vcf.gz
            rm temp_renamed.vcf.gz header.txt
        else
            mv temp_renamed.vcf.gz ~{prefix}.vcf.gz
        fi
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: 8 * ceil(size(vcf, "GB")) + 5,
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

task ResetVcfFilters {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x FILTER ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reset_vcf = "~{prefix}.vcf.gz"
        File reset_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task RevertSymbolicAlleles {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/helpers/revert_symbalts.py \
            ~{annotated_vcf} \
            ~{original_vcf} \
        | bcftools view -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.vcf.gz"
        File reverted_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotated_vcf, "GB") + size(original_vcf, "GB")) + 10,
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

task SplitQueryVcf {
    input {
        File vcf
        File? vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp.bed
        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed
        set +o pipefail
        head -1 ~{prefix}.bed > header
        set -o pipefail
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed )> ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed )> ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' ~{prefix}.bed )> ~{prefix}.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' ~{prefix}.bed )> ~{prefix}.INV_CPX.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' ~{prefix}.bed )> ~{prefix}.BND_CTX.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del_bed = "~{prefix}.DEL.bed"
        File dup_bed = "~{prefix}.DUP.bed"
        File ins_bed = "~{prefix}.INS.bed"
        File inv_bed = "~{prefix}.INV_CPX.bed"
        File bnd_bed = "~{prefix}.BND_CTX.bed"
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

task SplitVcfIntoShards {
  input {
    File input_vcf
    File input_vcf_idx
    Int variants_per_shard
    String prefix
    String docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    mkdir chunks
    bcftools view -h ~{input_vcf} > chunks/header.vcf
    bcftools view -H ~{input_vcf} | split -l ~{variants_per_shard} - chunks/body_

    for body in chunks/body_*; do
      chunk_name=chunks/~{prefix}_$(basename "$body")
      cat chunks/header.vcf "$body" | bgzip -c > "${chunk_name}.vcf.gz"
      tabix -p vcf -f "${chunk_name}.vcf.gz"
    done
 
   >>>

  output {
    Array[File] split_vcfs = glob("chunks/*.vcf.gz")
    Array[File] split_vcf_idxes = glob("chunks/*.vcf.gz.tbi")
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(size(input_vcf, "GB")) + 5,
    disk_gb: 2 * ceil(size(input_vcf,"GB")*2.5) + 10,
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

task StripGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -G ~{vcf} -Oz -o ~{prefix}.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File stripped_vcf = "~{prefix}.vcf.gz"
        File stripped_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task SubsetTsvToContig {
    input {
        File tsv
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -v contig="~{contig}" '$1 == contig' ~{tsv} > ~{prefix}.tsv
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
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

task SubsetVcfBySize {
    input {
        File vcf
        File vcf_idx
        String locus
        Int? min_size
        Int? max_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String size_filter = if defined(min_size) && defined(max_size) then "abs(INFO/SVLEN)>=~{min_size} && abs(INFO/SVLEN)<=~{max_size}" else if defined(min_size) then "abs(INFO/SVLEN)>=~{min_size}" else if defined(max_size) then "abs(INFO/SVLEN)<=~{max_size}" else "1==1"

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} --regions "~{locus}" --include "~{size_filter}" | bgzip > "~{prefix}.vcf.gz"
        
        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task SubsetVcfToCalled {
    input {
        File vcf
        File? vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if ~{!defined(vcf_idx)}; then
            tabix -p vcf ~{vcf}
        fi

        bcftools view \
            -i 'ALT != "."' \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task SubsetVcfToContig {
    input {
        File vcf
        File vcf_idx
        String contig
        String? args_string
        Boolean drop_genotypes = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -r ~{contig} \
            ~{if drop_genotypes then "-G" else ""} \
            ~{if defined(args_string) then args_string else ""} \
            ~{vcf} \
            -Oz -o ~{prefix}.~{contig}.vcf.gz
        tabix -p vcf -f ~{prefix}.~{contig}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.~{contig}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
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

task SubsetVcfToSampleList {
    input {
        File vcf
        File vcf_idx
        Array[String] samples
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > samples.txt <<EOF
~{sep='\n' samples}
EOF

        bcftools view \
            --samples-file samples.txt \
            --min-ac 1 \
            --regions ~{contig} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
