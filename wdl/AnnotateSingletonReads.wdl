version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow AnnotateSingletonReads {
    input {
        File vcf
        File vcf_index
        String pipeline_docker
        
        String output_suffix = "_singleton_filtered"
        Int? variants_per_shard

        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_postprocess
    }

    String prefix = basename(vcf, ".vcf.gz")

    Int variants_per_shard_eff = select_first([variants_per_shard, 1000000000])

    call Helpers.SplitVcfIntoShards {
        input:
            input_vcf = vcf,
            input_vcf_index = vcf_index,
            variants_per_shard = variants_per_shard_eff,
            output_prefix = prefix,
            docker_image = pipeline_docker
    }

    scatter (shard in zip(SplitVcfIntoShards.split_vcfs, SplitVcfIntoShards.split_vcf_indexes)) {
        String shard_prefix = basename(shard.left, ".vcf.gz")
        
        call FilterSingletonReads {
            input:
                vcf = shard.left,
                vcf_index = shard.right,
                prefix = shard_prefix,
                output_suffix = output_suffix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    call Helpers.ConcatVcfs as ConcatFilteredVcfs {
        input:
            vcfs = FilterSingletonReads.filtered_vcf,
            vcfs_idx = FilterSingletonReads.filtered_vcf_index,
            outfile_prefix = prefix,
            docker_image = pipeline_docker
    }

    call PostprocessVcf {
        input:
            vcf = ConcatFilteredVcfs.concat_vcf,
            vcf_index = ConcatFilteredVcfs.concat_vcf_idx,
            prefix = prefix,
            output_suffix = output_suffix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_postprocess
    }

    output {
        File singleton_filtered_vcf = PostprocessVcf.final_vcf
        File singleton_filtered_vcf_index = PostprocessVcf.final_vcf_index
    }
}

task FilterSingletonReads {
    input {
        File vcf
        File vcf_index
        String prefix
        String output_suffix
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

vcf_out = VariantFile("~{prefix}~{output_suffix}.vcf", 'w', header=header)

samples = list(vcf_in.header.samples)

for rec in vcf_in:
    is_suspicious = False
    if 'AC' in rec.info and rec.info['AC'][0] <= 2:  # Only check variants with AC <= 2
        alt_AD_in_called_samples = []
        for sample in samples:  # Add alt_AD for each sample if it's >0
            if 'AD' in rec.samples[sample]:
                ad_values = rec.samples[sample]['AD']
                if ad_values is not None and len(ad_values) == 2:
                    alt_AD = ad_values[1]
                    if alt_AD is not None and alt_AD > 0:
                        alt_AD_in_called_samples.append(alt_AD)
        
        if len(alt_AD_in_called_samples) == 1 and alt_AD_in_called_samples[0] == 1:
            is_suspicious = True
    
    if is_suspicious:
        rec.filter.add("SINGLE_READ_SUPPORT")
    
    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
EOF

        bgzip "~{prefix}~{output_suffix}.vcf"
        tabix "~{prefix}~{output_suffix}.vcf.gz"
    >>>

    output {
        File filtered_vcf = "~{prefix}~{output_suffix}.vcf.gz"
        File filtered_vcf_index = "~{prefix}~{output_suffix}.vcf.gz.tbi"
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

task PostprocessVcf {
    input {
        File vcf
        File vcf_index
        String prefix
        String output_suffix
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        # Simply copy the filtered VCF as final output
        cp ~{vcf} ~{prefix}~{output_suffix}.final.vcf.gz
        cp ~{vcf_index} ~{prefix}~{output_suffix}.final.vcf.gz.tbi
    >>>

    output {
        File final_vcf = "~{prefix}~{output_suffix}.final.vcf.gz"
        File final_vcf_index = "~{prefix}~{output_suffix}.final.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB")) + 5,
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
