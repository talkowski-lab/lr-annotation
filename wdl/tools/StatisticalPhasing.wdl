version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow StatisticalPhasing {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix

        Int operation
        String weight_tag
        Int is_weight_format_field
        Float default_weight
        Boolean do_shapeit5
        Boolean remove_duplicates_by_phased_fraction

        String variant_filter_args = "-i 'MAC>=2'"
        String filter_common_args = "-i 'MAF>=0.001'"
        String chunk_extra_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"
        String shapeit4_extra_args = "--thread $(nproc) --use-PS 0.0001"
        String shapeit5_extra_args =  "--thread $(nproc)"

        File genetic_maps_tsv
        File fix_variant_collisions_java

        String utils_docker
        String pysam_docker

        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_uqids
        RuntimeAttr? runtime_attr_remove_duplicates
        RuntimeAttr? runtime_attr_fill_vcf_tags
        RuntimeAttr? runtime_attr_fix_variant_collisions
        RuntimeAttr? runtime_attr_create_shapeit_chunks
        RuntimeAttr? runtime_attr_shapeit4_all
        RuntimeAttr? runtime_attr_filter_common
        RuntimeAttr? runtime_attr_shapeit4_common
        RuntimeAttr? runtime_attr_ligate_vcfs
        RuntimeAttr? runtime_attr_shapeit5_rare
        RuntimeAttr? runtime_attr_concat_shapeit5
        RuntimeAttr? runtime_attr_annotate_shapeit
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    call SplitAndFilterVcf {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = "~{prefix}.filtered",
            filter_args = variant_filter_args,
            runtime_attr_override = runtime_attr_filter_vcf
    }

    call EnsureUniqueIDsAndNormalize {
        input:
            vcf = SplitAndFilterVcf.filtered_vcf,
            vcf_idx = SplitAndFilterVcf.filtered_vcf_idx,
            prefix = "~{prefix}.filtered.uqids",
            docker = pysam_docker,
            runtime_attr_override = runtime_attr_uqids
    }

    if (remove_duplicates_by_phased_fraction) {
        call RemoveDuplicatesByPhasedFraction {
            input:
                vcf = EnsureUniqueIDsAndNormalize.uqids_vcf,
                vcf_idx = EnsureUniqueIDsAndNormalize.uqids_vcf_idx,
                prefix = "~{prefix}.filtered.uqids",
                docker = pysam_docker,
                runtime_attr_override = runtime_attr_remove_duplicates
        }
    }

    call FillVcfTags {
        input:
            vcf = select_first([RemoveDuplicatesByPhasedFraction.dups_removed_vcf, EnsureUniqueIDsAndNormalize.uqids_vcf]),
            vcf_idx = select_first([RemoveDuplicatesByPhasedFraction.dups_removed_vcf_idx, EnsureUniqueIDsAndNormalize.uqids_vcf_idx]),
            prefix = "~{prefix}.preprocessed",
            runtime_attr_override = runtime_attr_fill_vcf_tags
    }

    call FixVariantCollisions {
        input:
            phased_vcf = FillVcfTags.filled_tag_vcf,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            default_weight = default_weight,
            prefix = "~{prefix}.collisionless",
            runtime_attr_override = runtime_attr_fix_variant_collisions
    }

    call CreateShapeitChunks {
        input:
            vcf = FixVariantCollisions.phased_collisionless_vcf,
            vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
            region = contig,
            extra_args = chunk_extra_args,
            prefix = "~{prefix}.chunked",
            runtime_attr_override = runtime_attr_create_shapeit_chunks
    }

    Array[String] common_regions = read_lines(CreateShapeitChunks.common_chunks)

    scatter (i in range(length(common_regions))) {
        if (!do_shapeit5) {
            call Shapeit4 as Shapeit4All {
                input:
                    vcf = FixVariantCollisions.phased_collisionless_vcf,
                    vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
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
                    vcf = FixVariantCollisions.phased_collisionless_vcf,
                    vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
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
                    vcf = FixVariantCollisions.phased_collisionless_vcf,
                    vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
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

        call ConcatVcfs as ConcatShapeit5 {
            input:
                vcfs = flatten([Shapeit5Rare.phased_vcf]),
                vcf_idxs = flatten([Shapeit5Rare.phased_vcf_idx]),
                prefix = "~{prefix}.phased.concat",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shapeit5
        }
    }

    call AnnotateSHAPEITVcf {
        input:
            shapeit_input_vcf = FixVariantCollisions.phased_collisionless_vcf,
            shapeit_input_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
            shapeit_output_vcf = select_first([ConcatShapeit5.concat_vcf, LigateScaffold.ligated_vcf]),
            shapeit_output_vcf_idx = select_first([ConcatShapeit5.concat_vcf_idx, LigateScaffold.ligated_vcf_idx]),
            prefix = "~{prefix}.phased.annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_annotate_shapeit
    }

    output {
        File uqids_split_vcf = EnsureUniqueIDsAndNormalize.uqids_vcf
        File uqids_split_vcf_idx = EnsureUniqueIDsAndNormalize.uqids_vcf_idx
        File? removed_duplicates_split_vcf = select_first([FillVcfTags.filled_tag_vcf, RemoveDuplicatesByPhasedFraction.dups_removed_vcf])
        File? removed_duplicates_split_vcf_idx = select_first([FillVcfTags.filled_tag_vcf_idx, RemoveDuplicatesByPhasedFraction.dups_removed_vcf_idx])
        File collisionless_split_vcf = FixVariantCollisions.phased_collisionless_vcf
        File collisionless_split_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx
        File shapeit_phased_vcf = select_first([ConcatShapeit5.concat_vcf, LigateScaffold.ligated_vcf])
        File shapeit_phased_vcf_idx = select_first([ConcatShapeit5.concat_vcf_idx, LigateScaffold.ligated_vcf_idx])
        File shapeit_phased_annotated_vcf = AnnotateSHAPEITVcf.annotated_phased_vcf
        File shapeit_phased_annotated_vcf_idx = AnnotateSHAPEITVcf.annotated_phased_vcf_idx
    }
}

task SplitAndFilterVcf {
    input {
        File vcf
        File vcf_idx
        String prefix
        String filter_args
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm -m-any ~{vcf} | bcftools view ~{filter_args} -Oz -o ~{prefix}.vcf.gz
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task EnsureUniqueIDsAndNormalize {
  input {
      File vcf
      File vcf_idx
      String prefix
      String docker
      RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    python3 <<CODE

    #!/usr/bin/env python3
    import pysam

    def remove_longest_shared_suffix(s1, s2):
        i = 0
        l1 = len(s1)
        l2 = len(s2)
        while i < l1 - 1 and i < l2 - 1:
            if s1[(l1 - i - 1)] == s2[(l2 - i - 1)]:
                i += 1
            else:
                break
        return s1[:(l1 - i)], s2[:(l2 - i)]

    def remove_longest_shared_prefix(s1, s2):
        i = 0
        l1 = len(s1)
        l2 = len(s2)
        while i < l1 - 1 and i < l2 - 1:
            if s1[i] == s2[i]:
                i += 1
            else:
                break
        return s1[i:], s2[i:], i

    def ensure_unique_ids_and_normalize(input_vcf, output_vcf):
        vcf_in = pysam.VariantFile(input_vcf, "r")
        vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)
        id_dict = dict()
        for record in vcf_in:
            if record.id not in id_dict:
                id_dict[record.id] = 1
            else:
                id_dict[record.id] = id_dict[record.id] + 1
            assert len(record.alts) == 1
            new_ref, new_alt, i = remove_longest_shared_prefix(*remove_longest_shared_suffix(record.ref, record.alts[0]))
            new_record = record.copy()
            new_record.id = record.id + "_" + str(id_dict[record.id])
            new_record.pos = record.pos + i
            new_record.ref = new_ref
            new_record.alts = [new_alt]
            vcf_out.write(new_record)
        vcf_in.close()
        vcf_out.close()
        print(f"✅ Unique ID VCF written to: {output_vcf}")

    input_vcf = "~{vcf}"
    output_vcf = "unsorted.vcf"
    ensure_unique_ids_and_normalize(input_vcf, output_vcf)

    CODE

    bcftools sort "unsorted.vcf" -Oz -o "~{prefix}.vcf.gz"
    bcftools index -t "~{prefix}.vcf.gz"

  >>>

  output {
    File uqids_vcf = "~{prefix}.vcf.gz"
    File uqids_vcf_idx =  "~{prefix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf, "GiB")*3),
    disk_gb: 15 + ceil(size(vcf, "GiB")*3),
    boot_disk_gb: 10,
    preemptible_tries: 1,
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

task RemoveDuplicatesByPhasedFraction {
  input {
      File vcf
      File vcf_idx
      String prefix
      String docker
      RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    python3 <<CODE

    #!/usr/bin/env python3
    import pysam

    def annotate_phased_fraction(input_vcf, output_vcf):
        vcf_in = pysam.VariantFile(input_vcf, "r")
        vcf_in.header.add_meta("INFO", items=[("ID", "PHYS_PHASE_FRAC"),
                                              ("Number", "1"),
                                              ("Type", 'Float'),
                                              ("Description", "Fraction of genotyped samples physically phased")])
        vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)
        v_to_max_phased_frac_and_id_dict = dict()
        remove_ids = set()
        for record in vcf_in:
            v = (record.contig, record.pos, record.ref.upper(), tuple([a.upper() for a in record.alts]))
            phased_count = sum([s.phased and s["GT"] is not None and not None in s["GT"] for s in record.samples.values()])
            genotyped_count = sum([s["GT"] is not None and not None in s["GT"] for s in record.samples.values()])
            if genotyped_count == 0:
                phased_frac = 0
            else:
                phased_frac = round(phased_count / genotyped_count, 4)
            if v not in v_to_max_phased_frac_and_id_dict:
                v_to_max_phased_frac_and_id_dict[v] = (phased_frac, record.id)
            elif phased_frac < v_to_max_phased_frac_and_id_dict[v][0]:
                remove_ids.add(record.id)
            elif phased_frac > v_to_max_phased_frac_and_id_dict[v][0]:
                remove_ids.add(v_to_max_phased_frac_and_id_dict[v][1])
                v_to_max_phased_frac_and_id_dict[v] = (phased_frac, record.id)
            elif phased_frac == v_to_max_phased_frac_and_id_dict[v][0]:
                remove_ids.add(record.id) # just keep the first variant
            new_record = record.copy()
            new_record.info["PHYS_PHASE_FRAC"] = phased_frac
            new_record.ref = new_record.ref.upper()
            new_record.alts = [a.upper() for a in new_record.alts]
            vcf_out.write(new_record)
        vcf_in.close()
        vcf_out.close()
        print(f"✅ Annotated VCF written to: {output_vcf}")
        print("Number of duplicate variants to remove: " + str(len(remove_ids)))
        return remove_ids

    def remove_duplicate_variants(input_vcf, output_vcf, remove_ids):
        vcf_in = pysam.VariantFile(input_vcf, "r")
        vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)
        for record in vcf_in:
            if record.id not in remove_ids:
                vcf_out.write(record)
        vcf_in.close()
        vcf_out.close()
        print(f"✅ VCF with duplicate variants removed written to: {output_vcf}")

    input_vcf = "~{vcf}"
    intermediate_vcf = "~{prefix}.annotated.vcf"
    intermediate_vcf_gz = "~{prefix}.annotated.vcf.gz"
    output_vcf = "~{prefix}.annotated.dups_removed.vcf"

    remove_ids = annotate_phased_fraction(input_vcf, intermediate_vcf)
    pysam.tabix_compress(intermediate_vcf, intermediate_vcf_gz)
    pysam.tabix_index(intermediate_vcf_gz, preset="vcf")

    remove_duplicate_variants(intermediate_vcf_gz, output_vcf, remove_ids)

    CODE

    bgzip "~{prefix}.annotated.dups_removed.vcf"
    bcftools index -t "~{prefix}.annotated.dups_removed.vcf.gz"

  >>>

  output {
    File annotated_vcf = "~{prefix}.annotated.vcf.gz"
    File annotated_vcf_idx = "~{prefix}.annotated.vcf.gz.tbi"
    File dups_removed_vcf = "~{prefix}.annotated.dups_removed.vcf.gz"
    File dups_removed_vcf_idx =  "~{prefix}.annotated.dups_removed.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf, "GiB")*3),
    disk_gb: 15 + ceil(size(vcf, "GiB")*3),
    boot_disk_gb: 10,
    preemptible_tries: 1,
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

task FillVcfTags {
    input {
        File vcf
        File vcf_idx
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools +setGT ~{vcf} -- -t . -n 0p | \
            bcftools +fill-tags -- -t AF,AC,AN,MAF | \
            bcftools view -Oz -o ~{prefix}.vcf.gz

        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filled_tag_vcf = "~{prefix}.vcf.gz"
        File filled_tag_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        bcftools +fill-tags \
            -r ~{region} \
            ~{vcf} \
            -- \
            -t AF,AC,AN \
        | bcftools view \
            ~{filter_common_args} \
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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

        shapeit4.2 \
            --input ~{vcf} \
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "hangsuunc/shapeit5:v1"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit5:develop"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Boolean merge_sort = true
        String prefix = "concat"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCFS_FILE="~{write_lines(vcfs)}"

        bcftools concat \
            ~{if merge_sort then "--allow-overlaps" else ""} \
            --file-list ${VCFS_FILE} \
            -Oz -o "~{prefix}.vcf.gz"

        tabix -p vcf -f "~{prefix}.vcf.gz"
    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
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

task AnnotateSHAPEITVcf {
    input {
        File shapeit_input_vcf
        File shapeit_input_vcf_idx
        File shapeit_output_vcf
        File shapeit_output_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bcftools annotate -a ~{shapeit_output_vcf} ~{shapeit_input_vcf} -c FORMAT/GT,INFO -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

      >>>

    output {
        File annotated_phased_vcf = "~{prefix}.vcf.gz"
        File annotated_phased_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(shapeit_input_vcf, "GiB")*3) + ceil(size(shapeit_output_vcf, "GiB")*3),
        disk_gb: 15 + ceil(size(shapeit_input_vcf, "GiB")*3) + ceil(size(shapeit_output_vcf, "GiB")*3),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
