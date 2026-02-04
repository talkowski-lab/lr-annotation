version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow SubsetTRGTToCatalog {
    input {
        File trgt_merged_vcf
        File trgt_merged_vcf_idx
        File trgt_catalog_bed_gz
        Array[String] contigs
        String prefix
        
        String utils_docker
        
        RuntimeAttr? runtime_attr_extract_ids
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_filter_to_catalog
        RuntimeAttr? runtime_attr_concat_vcf
    }
    
    call ExtractCatalogIDs {
        input:
            catalog_bed_gz = trgt_catalog_bed_gz,
            contigs = contigs,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_extract_ids
    }
    
    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = trgt_merged_vcf,
                vcf_idx = trgt_merged_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
        
        call FilterToCatalog {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                catalog_ids = ExtractCatalogIDs.catalog_ids,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_to_catalog
        }
    }
    
    call Helpers.ConcatVcfs {
        input:
            vcfs = FilterToCatalog.filtered_vcf,
            vcfs_idx = FilterToCatalog.filtered_vcf_idx,
            prefix = prefix + ".catalog_subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }
    
    output {
        File catalog_subset_vcf = ConcatVcfs.concat_vcf
        File catalog_subset_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task ExtractCatalogIDs {
    input {
        File catalog_bed_gz
        Array[String] contigs
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -euo pipefail
        
        cat > contigs.txt <<EOF
~{sep="\n" contigs}
EOF
        
        zcat ~{catalog_bed_gz} | \
            awk -F'\t' 'BEGIN {
                while ((getline < "contigs.txt") > 0) {
                    contig_filter[$0] = 1;
                }
            }
            $1 in contig_filter {
                split($4, fields, ";");
                for (i in fields) {
                    if (fields[i] ~ /^ID=/) {
                        sub(/^ID=/, "", fields[i]);
                        print fields[i];
                        break;
                    }
                }
            }' | \
            sort -u > catalog_ids.txt
    >>>
    
    output {
        File catalog_ids = "catalog_ids.txt"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(catalog_bed_gz, "GB")) + 10,
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

task FilterToCatalog {
    input {
        File vcf
        File vcf_idx
        File catalog_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -euo pipefail
        
        cp ~{catalog_ids} catalog_ids.txt
        
        cat > filter_script.awk <<'EOF'
BEGIN {
    while ((getline < "catalog_ids.txt") > 0) {
        allowed_ids[$0] = 1;
    }
}
!/^#/ {
    trid = "";
    n = split($8, info_fields, ";");
    for (i = 1; i <= n; i++) {
        if (info_fields[i] ~ /^TRID=/) {
            sub(/^TRID=/, "", info_fields[i]);
            trid = info_fields[i];
            break;
        }
    }
    if (trid != "" && trid in allowed_ids) {
        print;
    }
}
EOF
        
        bcftools view -h ~{vcf} > header.vcf
        bcftools view -H ~{vcf} | \
            awk -f filter_script.awk | \
            cat header.vcf - | \
            bgzip -c > ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>
    
    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
