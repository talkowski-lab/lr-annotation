version 1.0

import "general/Structs.wdl"
import "general/Helpers.wdl" as Helpers

workflow AnnotateSingletonReads {
    input {
        File vcf
        File vcf_index
        String pipeline_docker
        String prefix
        Array[String] contigs
        
        Int? variants_per_shard

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_split
        RuntimeAttr? runtime_attr_populate_tags
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat
    }

    Int variants_per_shard_eff = select_first([variants_per_shard, 1000000000])

    call SubsetVcfToContigs {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            contigs = contigs,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_subset
    }

    call Helpers.SplitVcfIntoShards {
        input:
            input_vcf = SubsetVcfToContigs.subset_vcf,
            input_vcf_index = SubsetVcfToContigs.subset_vcf_index,
            variants_per_shard = variants_per_shard_eff,
            output_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_split
    }

    scatter (shard in zip(SplitVcfIntoShards.split_vcfs, SplitVcfIntoShards.split_vcf_indexes)) {
        String shard_prefix = basename(shard.left, ".vcf.gz")
        
        call PopulateTags {
            input:
                vcf = shard.left,
                vcf_index = shard.right,
                prefix = shard_prefix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_populate_tags
        }
        
        call FilterSingletonReads {
            input:
                vcf = PopulateTags.tagged_vcf,
                vcf_index = PopulateTags.tagged_vcf_index,
                prefix = shard_prefix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    call Helpers.ConcatVcfs as ConcatFilteredVcfs {
        input:
            vcfs = FilterSingletonReads.filtered_vcf,
            vcfs_idx = FilterSingletonReads.filtered_vcf_index,
            outfile_prefix = prefix,
            docker_image = pipeline_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File singleton_filtered_vcf = ConcatFilteredVcfs.concat_vcf
        File singleton_filtered_vcf_index = ConcatFilteredVcfs.concat_vcf_idx
    }
}

task PopulateTags {
    input {
        File vcf
        File vcf_index
        String prefix
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        bcftools +fill-tags ~{vcf} -Oz -o ~{prefix}.tagged.vcf.gz -- -t AC
        tabix -p vcf ~{prefix}.tagged.vcf.gz
    >>>

    output {
        File tagged_vcf = "~{prefix}.tagged.vcf.gz"
        File tagged_vcf_index = "~{prefix}.tagged.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB") * 2) + 10,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FilterSingletonReads {
    input {
        File vcf
        File vcf_index
        String prefix
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        python3 <<'EOF'
import sys
from pysam import VariantFile

vcf_in = VariantFile("~{vcf}")

header = vcf_in.header
header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by only one read in a single sample")

vcf_out = VariantFile("~{prefix}.filtered.vcf", 'w', header=header)

samples = list(vcf_in.header.samples)

for rec in vcf_in:
    is_suspicious = False
    
    if rec.info['AC'][0] <= 2:
        alt_AD_in_called_samples = []
        for sample in samples:
            if 'AD' in rec.samples[sample]:
                if len(rec.samples[sample]['AD']) == 2:
                    alt_AD = rec.samples[sample]['AD'][1]
                    if alt_AD > 0:
                        alt_AD_in_called_samples.append(alt_AD)
        if len(alt_AD_in_called_samples)==1 and alt_AD_in_called_samples[0]==1:
            is_suspicious = True
    
    if is_suspicious:
        rec.filter.add("SINGLE_READ_SUPPORT")
    
    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
EOF

        bgzip "~{prefix}.filtered.vcf"
        tabix "~{prefix}.filtered.vcf.gz"
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 5,
        disk_gb: ceil(size(vcf, "GB")) + 20,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToContigs {
    input {
        File vcf
        File vcf_index
        Array[String] contigs
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} --regions ~{sep=',' contigs} -Oz -o ~{prefix}.subset.vcf.gz
        tabix -p vcf ~{prefix}.subset.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.subset.vcf.gz"
        File subset_vcf_index = "~{prefix}.subset.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB")) + 20,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
