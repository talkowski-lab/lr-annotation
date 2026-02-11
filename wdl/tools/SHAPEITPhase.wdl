version 1.0

import "../utils/Structs.wdl" as Structs
import "../utils/Helpers.wdl" as Helpers

workflow SHAPEITPhase {
    input {
        File joint_vcf
        File joint_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        File genetic_maps_tsv
        String chromosome
        String region
        String prefix

        Int shard_size = 1000000
        String variant_filter_args = "-i 'MAC>=2'"

        File fix_variant_collisions_java
        Int operation = 1
        String weight_tag = "SCORE"
        Int is_weight_format_field = 0
        Float default_weight = 0.5

        String chunk_extra_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"

        Boolean do_shapeit5 = true
        String shapeit4_extra_args = "--thread $(nproc) --use-PS 0.0001"
        String shapeit5_extra_args =  "--thread $(nproc)"
        String filter_common_args = "-i 'MAF>=0.001'"

        String docker = "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_fix_variant_collisions
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_create_shapeit_chunks
        RuntimeAttr? runtime_attr_shapeit4
        RuntimeAttr? runtime_attr_filter_common
        RuntimeAttr? runtime_attr_ligate_vcfs
        RuntimeAttr? runtime_attr_shapeit5_rare
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    call SubsetVcfToRegion {
        input:
            vcf = joint_vcf,
            vcf_idx = joint_vcf_idx,
            region = region,
            prefix = "~{prefix}.region",
            runtime_attr_override = runtime_attr_subset_vcf
    }

    scatter (s in range(length(SubsetVcfToRegion.shard_regions))) {
        String shard_region = SubsetVcfToRegion.shard_regions[s]

        call SubsetVcfToRegion as SubsetShard {
            input:
                vcf = SubsetVcfToRegion.subset_vcf,
                vcf_idx = SubsetVcfToRegion.subset_vcf_idx,
                region = shard_region,
                prefix = "~{prefix}.shard-~{s}",
                shard_size = shard_size,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call SplitAndFilterVcf {
            input:
                vcf = SubsetShard.subset_vcf,
                vcf_idx = SubsetShard.subset_vcf_idx,
                prefix = "~{prefix}.filtered.shard-~{s}",
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                filter_args = variant_filter_args,
                runtime_attr_override = runtime_attr_filter_vcf
        }

        call FixVariantCollisions {
            input:
                phased_vcf = SplitAndFilterVcf.filtered_vcf,
                fix_variant_collisions_java = fix_variant_collisions_java,
                operation = operation,
                weight_tag = weight_tag,
                is_weight_format_field = is_weight_format_field,
                default_weight = default_weight,
                prefix = "~{prefix}.collisionless.shard-~{s}",
                runtime_attr_override = runtime_attr_fix_variant_collisions
        }
    }

    call Helpers.ConcatVcfs as ConcatShards {
        input:
            vcfs = FixVariantCollisions.phased_collisionless_vcf,
            vcf_idxs = FixVariantCollisions.phased_collisionless_vcf_idx,
            merge_sort = false,
            prefix = "~{prefix}.prepared",
            docker = docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    call CreateShapeitChunks {
        input:
            vcf = ConcatShards.concat_vcf,
            vcf_idx = ConcatShards.concat_vcf_idx,
            region = region,
            extra_args = chunk_extra_args,
            runtime_attr_override = runtime_attr_create_shapeit_chunks
    }

    Array[String] common_regions = read_lines(CreateShapeitChunks.common_chunks)

    scatter (i in range(length(common_regions))) {
        if (!do_shapeit5) {
            call Shapeit4 as Shapeit4All {
                input:
                    vcf = ConcatShards.concat_vcf,
                    vcf_idx = ConcatShards.concat_vcf_idx,
                    genetic_map = genetic_maps_dict[chromosome],
                    region = common_regions[i],
                    prefix = "~{prefix}.phased",
                    extra_args = shapeit4_extra_args,
                    runtime_attr_override = runtime_attr_shapeit4
            }
        }
        if (do_shapeit5) {
            call FilterCommon {
                input:
                    vcf = ConcatShards.concat_vcf,
                    vcf_idx = ConcatShards.concat_vcf_idx,
                    prefix = "~{prefix}.common",
                    region = common_regions[i],
                    filter_common_args = filter_common_args,
                    runtime_attr_override = runtime_attr_filter_common
            }

            call Shapeit4 as Shapeit4Common {
                input:
                    vcf = FilterCommon.common_vcf,
                    vcf_idx = FilterCommon.common_vcf_idx,
                    genetic_map = genetic_maps_dict[chromosome],
                    region = common_regions[i],
                    prefix = "~{prefix}.phased",
                    extra_args = shapeit4_extra_args,
                    runtime_attr_override = runtime_attr_shapeit4
            }
        }
    }

    call LigateVcfs as LigateScaffold {
        input:
            vcfs = select_all(flatten([Shapeit4All.phased_vcf, Shapeit4Common.phased_vcf])),
            vcf_idxs = select_all(flatten([Shapeit4All.phased_vcf_idx, Shapeit4Common.phased_vcf_idx])),
            prefix = "~{prefix}.phased.ligated",
            runtime_attr_override = runtime_attr_ligate_vcfs
    }

    if (do_shapeit5) {
        Array[String] rare_regions = read_lines(CreateShapeitChunks.rare_chunks)
        scatter (i in range(length(rare_regions))) {
            call Shapeit5Rare {
                input:
                    vcf = ConcatShards.concat_vcf,
                    vcf_idx = ConcatShards.concat_vcf_idx,
                    scaffold_vcf = LigateScaffold.ligated_vcf,
                    scaffold_vcf_idx = LigateScaffold.ligated_vcf_idx,
                    genetic_map = genetic_maps_dict[chromosome],
                    region = rare_regions[i],
                    scaffold_region = common_regions[i],
                    prefix = "~{prefix}.phased.chunk-~{i}",
                    extra_args = shapeit5_extra_args,
                    runtime_attr_override = runtime_attr_shapeit5_rare
            }
        }

        call Helpers.ConcatVcfs as ConcatShapeit5 {
            input:
                vcfs = flatten([Shapeit5Rare.phased_vcf]),
                vcf_idxs = flatten([Shapeit5Rare.phased_vcf_idx]),
                merge_sort = false,
                prefix = "~{prefix}.phased.concat",
                docker = docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        File shapeit_phased_vcf = select_first([ConcatShapeit5.concat_vcf, LigateScaffold.ligated_vcf])
        File shapeit_phased_vcf_idx = select_first([ConcatShapeit5.concat_vcf_idx, LigateScaffold.ligated_vcf_idx])
    }
}

task SubsetVcfToRegion {
    input {
        File vcf
        File vcf_idx
        String region
        String prefix
        Int shard_size = 1000000
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([vcf, vcf_idx], "GB")) + 20

    command <<<
        set -euxo pipefail

        bcftools view --no-version ~{vcf} --regions ~{region} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

        python3 <<CODE
region = "~{region}"
chrom, span = region.split(":")
start, end = map(int, span.split("-"))
bin_size = ~{shard_size}

intervals = []
pos = start
while pos < end:
    interval_end = min(pos + bin_size, end)
    intervals.append(f"{chrom}:{pos}-{interval_end}")
    pos = interval_end

with open("shards.txt", "w") as f:
    for interval in intervals:
        f.write(interval + "\n")
CODE
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
        Array[String] shard_regions = read_lines("shards.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task SplitAndFilterVcf {
    input {
        File vcf
        File vcf_idx
        String prefix
        File reference_fasta
        File reference_fasta_fai
        String filter_args
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 4 * ceil(size(vcf, "GB")) + 10

    command <<<
        set -euxo pipefail

        bcftools norm --no-version -m-any -N -f ~{reference_fasta} ~{vcf} | \
            bcftools +setGT --no-version -- -t . -n 0p | \
            bcftools +fill-tags --no-version -- -t AF,AC,AN,MAF | \
            bcftools view --no-version ~{filter_args} -Oz -o ~{prefix}.vcf.gz
        
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task FixVariantCollisions {
    input {
        File phased_vcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "SCORE"         # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight = 0.5          # default weight if the weight field is not found
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 100 + 4 * (ceil(size(phased_vcf, "GiB")))

    command <<<
        set -euxo pipefail

        java ~{fix_variant_collisions_java} \
            ~{phased_vcf} \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools +setGT --no-version collisionless.vcf -- -t . -n 0p | \
            bcftools +fill-tags --no-version -Oz -o ~{prefix}.phased.collisionless.vcf.gz -- -t AF,AC,AN
        # use vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{prefix}.phased.collisionless.vcf.gz
    >>>

    output {
        File phased_collisionless_vcf = "~{prefix}.phased.collisionless.vcf.gz"
        File phased_collisionless_vcf_idx = "~{prefix}.phased.collisionless.vcf.gz.tbi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"     # needs Java + bcftools
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task CreateShapeitChunks {
    input {
        File vcf
        File vcf_idx
        String region
        String extra_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([vcf, vcf_idx], "GiB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > common.chunks.regions.txt
        cut -f 4 chunks.txt > rare.chunks.regions.txt
    >>>

    output {
        File common_chunks = "common.chunks.regions.txt"
        File rare_chunks = "rare.chunks.regions.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task BcftoolsConcatNaive {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        bcftools concat --no-version ~{sep=" " vcfs} --naive -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File concatenated_vcf = "~{prefix}.vcf.gz"
        File concatenated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task FilterCommon {
    input {
        File vcf
        File vcf_idx
        String prefix
        String region
        String filter_common_args
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * ceil(size(vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # filter to common
        bcftools +fill-tags --no-version -r ~{region} ~{vcf} -- -t AF,AC,AN | \
            bcftools view ~{filter_common_args} \
                -Ob -o ~{prefix}.common.bcf
        bcftools index ~{prefix}.common.bcf
    >>>

    output {
        File common_vcf = "~{prefix}.common.bcf"
        File common_vcf_idx = "~{prefix}.common.bcf.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task Shapeit4 {
    input {
        File vcf
        File vcf_idx
        File genetic_map
        String region
        String prefix
        String extra_args = "--thread $(nproc) --use-PS 0.0001"
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB"))

    command <<<
        set -euxo pipefail

        shapeit4.2 --input ~{vcf} \
                --map ~{genetic_map} \
                --region ~{region} \
                --sequencing \
                --output ~{prefix}.bcf \
                ~{extra_args}
        bcftools index ~{prefix}.bcf
    >>>

    output{
        File phased_vcf = "~{prefix}.bcf"
        File phased_vcf_idx = "~{prefix}.bcf.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 16,
        mem_gb: 16,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task LigateVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.bcf
        bcftools +fill-tags --no-version ~{prefix}.bcf \
            -Oz -o ~{prefix}.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf = "~{prefix}.vcf.gz"
        File ligated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task Shapeit5Rare {
    input {
        File vcf
        File vcf_idx
        File scaffold_vcf
        File scaffold_vcf_idx
        File genetic_map
        String region
        String scaffold_region
        String prefix
        String extra_args = "--thread $(nproc)"
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB")) + ceil(size(scaffold_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # we only need to fill rare in input (common in scaffold should have been imputed or filled previously);
        # this also fills common in input, but those records will be ignored by Shapeit5
        bcftools +setGT --no-version ~{vcf} -- -t . -n 0p | \
            bcftools +fill-tags --no-version -Ob -o input.bcf -- -t AF,AC,AN
        bcftools index input.bcf

        /shapeit5/phase_rare --input input.bcf \
            --scaffold ~{scaffold_vcf} \
            --map ~{genetic_map} \
            --input-region ~{region} \
            --scaffold-region ~{scaffold_region} \
            --output phased.bcf \
            ~{extra_args}

        bcftools +fill-tags phased.bcf --no-version -Oz -o ~{prefix}.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{prefix}.vcf.gz

    >>>

    output{
        File phased_vcf = "~{prefix}.vcf.gz"
        File phased_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 16,
        mem_gb: 32,
        disk_gb: disk_gb,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit5:develop"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}
