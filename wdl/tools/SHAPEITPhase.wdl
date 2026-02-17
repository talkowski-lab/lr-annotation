version 1.0

import "../utils/Structs.wdl" as Structs
import "../utils/Helpers.wdl" as Helpers

workflow SHAPEITPhase {
    input {
        File vcf
        File vcf_idx
        String contig
        String region
        String prefix

        Int operation
        String weight_tag
        Int is_weight_format_field
        Float default_weight
        Boolean do_shapeit5

        String variant_filter_args = "-i 'MAC>=2'"
        String filter_common_args = "-i 'MAF>=0.001'"
        String chunk_extra_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"
        String shapeit4_extra_args = "--thread $(nproc) --use-PS 0.0001"
        String shapeit5_extra_args =  "--thread $(nproc)"

        File ref_fa
        File ref_fai
        File genetic_maps_tsv
        File fix_variant_collisions_java

        Int filter_and_concat_shard_size
        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_fix_variant_collisions
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_create_shapeit_chunks
        RuntimeAttr? runtime_attr_shapeit4_all
        RuntimeAttr? runtime_attr_filter_common
        RuntimeAttr? runtime_attr_shapeit4_common
        RuntimeAttr? runtime_attr_ligate_vcfs
        RuntimeAttr? runtime_attr_shapeit5_rare
        RuntimeAttr? runtime_attr_concat_shapeit5
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    call CreateShards {
        input:
            region = region,
            bin_size = filter_and_concat_shard_size,
            pad_size = 0,
            prefix = "~{prefix}.shards",
            runtime_attr_override = runtime_attr_create_shards
    }

    scatter (i in range(length(CreateShards.shard_regions))) {
        String shard_region = CreateShards.shard_regions[i]

        call SubsetVcfToRegion as SubsetShard {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                region = shard_region,
                prefix = "~{prefix}.~{i}",
                runtime_attr_override = runtime_attr_subset_shard
        }

        call SplitAndFilterVcf {
            input:
                vcf = SubsetShard.subset_vcf,
                vcf_idx = SubsetShard.subset_vcf_idx,
                prefix = "~{prefix}.~{i}.filtered",
                ref_fa = ref_fa,
                ref_fai = ref_fai,
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
                prefix = "~{prefix}.~{i}.collisionless",
                runtime_attr_override = runtime_attr_fix_variant_collisions
        }
    }

    call Helpers.ConcatVcfs as ConcatShards {
        input:
            vcfs = FixVariantCollisions.phased_collisionless_vcf,
            vcf_idxs = FixVariantCollisions.phased_collisionless_vcf_idx,
            merge_sort = false,
            prefix = "~{prefix}.prepared",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_shards
    }

    call CreateShapeitChunks {
        input:
            vcf = ConcatShards.concat_vcf,
            vcf_idx = ConcatShards.concat_vcf_idx,
            region = region,
            extra_args = chunk_extra_args,
            prefix = "~{prefix}.chunked",
            runtime_attr_override = runtime_attr_create_shapeit_chunks
    }

    Array[String] common_regions = read_lines(CreateShapeitChunks.common_chunks)

    scatter (i in range(length(common_regions))) {
        if (!do_shapeit5) {
            call Shapeit4 as Shapeit4All {
                input:
                    vcf = ConcatShards.concat_vcf,
                    vcf_idx = ConcatShards.concat_vcf_idx,
                    genetic_map = genetic_maps_dict[contig],
                    region = common_regions[i],
                    prefix = "~{prefix}.phased",
                    extra_args = shapeit4_extra_args,
                    runtime_attr_override = runtime_attr_shapeit4_all
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
                    genetic_map = genetic_maps_dict[contig],
                    region = common_regions[i],
                    prefix = "~{prefix}.phased",
                    extra_args = shapeit4_extra_args,
                    runtime_attr_override = runtime_attr_shapeit4_common
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
                    genetic_map = genetic_maps_dict[contig],
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
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shapeit5
        }
    }

    output {
        File shapeit_phased_vcf = select_first([ConcatShapeit5.concat_vcf, LigateScaffold.ligated_vcf])
        File shapeit_phased_vcf_idx = select_first([ConcatShapeit5.concat_vcf_idx, LigateScaffold.ligated_vcf_idx])
    }
}

task CreateShards {
    input {
        String region
        Int bin_size
        Int pad_size
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python - --region ~{region} \
                 --bin_size ~{bin_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{prefix} \
                 <<-'EOF'
        import argparse

        def split_region(region):
            chromosome, span = region.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_region_to_intervals(region, bin_size, pad_size):
            chromo, start, end = split_region(region)
            bin_num = (end - start)//bin_size
            intervals = [(chromo, start, start + bin_size + pad_size)]
            for i in range(1, bin_num):
                start_pos = start + i*bin_size - pad_size
                end_pos = start_pos + bin_size + pad_size
                intervals.append((chromo, start_pos, end_pos))
            if end > start + bin_num*bin_size:
                intervals.append((chromo, start + bin_num*bin_size - pad_size, end))
            return(intervals)

        def write_output_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    l = "%s:%d-%d" % (item[0], item[1], item[2])
                    f.write(l+ "\n")

        def main():
            parser = argparse.ArgumentParser()
            parser.add_argument('--region', type=str)
            parser.add_argument('--output_file', type=str)
            parser.add_argument('--bin_size', type=int)
            parser.add_argument('--pad_size', type=int)
            args = parser.parse_args()

            intervals = split_region_to_intervals(args.region, args.bin_size, args.pad_size)
            write_output_file(intervals, args.output_file + ".txt")

        if __name__ == "__main__":
            main()
        EOF
    >>>

    output {
        Array[String] shard_regions = read_lines("~{prefix}.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task SubsetVcfToRegion {
    input {
        File vcf
        File vcf_idx
        String region
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -r ~{region} \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}
        
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
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
        File ref_fa
        File ref_fai
        String filter_args
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm -m-any -N -f ~{ref_fa} ~{vcf} | \
            bcftools +setGT -- -t . -n 0p | \
            bcftools +fill-tags -- -t AF,AC,AN,MAF | \
            bcftools view ~{filter_args} -Oz -o ~{prefix}.vcf.gz
        
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
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
        File phased_vcf
        File fix_variant_collisions_java
        Int operation                       # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag                   # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field          # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight                # default weight if the weight field is not found
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        java ~{fix_variant_collisions_java} \
            ~{phased_vcf} \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools +setGT collisionless.vcf -- -t . -n 0p | \
            bcftools +fill-tags -Oz -o ~{prefix}.phased.collisionless.vcf.gz -- -t AF,AC,AN
        
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
        disk_gb: 5 * ceil(size(phased_vcf, "GiB")) + 100,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"
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
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_args} \
            -O chunks.txt

        cut -f3 chunks.txt > ~{prefix}.common.txt
        cut -f4 chunks.txt > ~{prefix}.rare.txt
    >>>

    output {
        File common_chunks = "~{prefix}.common.txt"
        File rare_chunks = "~{prefix}.rare.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GiB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
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
        set -euo pipefail

        # filter to common
        bcftools +fill-tags -r ~{region} ~{vcf} -- -t AF,AC,AN | \
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
        max_retries: 0,
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
        set -euo pipefail

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
        max_retries: 0,
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
        set -euo pipefail

        ligate_static \
            --input ~{write_lines(vcfs)} \
            --output ~{prefix}.bcf
        
        bcftools +fill-tags \
            ~{prefix}.bcf \
            -Oz -o ~{prefix}.vcf.gz \
            -- -t AF,AC,AN
        
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
        max_retries: 0,
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
        set -euo pipefail

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
        max_retries: 0,
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
