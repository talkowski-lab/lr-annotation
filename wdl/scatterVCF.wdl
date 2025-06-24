version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow ScatterVCF {
    input {
        File file
        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/main/scripts/split_vcf_hail.py"
        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        Boolean localize_vcf
        Boolean split_by_chromosome
        Boolean split_into_shards 
        Boolean get_chromosome_sizes
        Boolean has_index=false
        Int n_shards=0
        Int records_per_shard=0
        String genome_build='GRCh38'
        RuntimeAttr? runtime_attr_split_by_chr
        RuntimeAttr? runtime_attr_split_into_shards
    }
    
    # shard the VCF (if not already sharded)
    if (split_by_chromosome) {
        if (!localize_vcf) {
            String vcf_uri = file
            if (get_chromosome_sizes) {
                call GetChromosomeSizes {
                    input:
                        vcf_file=vcf_uri,
                        has_index=select_first([has_index]),
                        sv_base_mini_docker=sv_base_mini_docker
                }
            }
        }

        Map[String, Array[String]] chromosomes_dict = {
            'GRCh38': ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
            'GRCh37': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']
        }
        Array[String] chromosomes = chromosomes_dict[genome_build]
        scatter (chromosome in chromosomes) {
            if (localize_vcf) {
                call SplitByChromosome {
                    input:
                        vcf_file=file,
                        chromosome=chromosome,
                        sv_base_mini_docker=sv_base_mini_docker,
                        runtime_attr_override=runtime_attr_split_by_chr
                }
            }
            if (!localize_vcf) {
                String vcf_uri = file
                Float input_size_ = if (get_chromosome_sizes) then 
                # chrom_length * ceil(n_samples*0.001) / 1000000
                    select_first([GetChromosomeSizes.contig_lengths])[chromosome] * ceil(select_first([GetChromosomeSizes.n_samples])*0.001) / 1000000 
                    else size(vcf_uri, 'GB')
                call SplitByChromosomeRemote {
                    input:
                        vcf_file=vcf_uri,
                        chromosome=chromosome,
                        input_size=input_size_,
                        has_index=select_first([has_index]),
                        sv_base_mini_docker=sv_base_mini_docker,
                        runtime_attr_override=runtime_attr_split_by_chr
                }
            }
            File splitChromosomeShards = select_first([SplitByChromosome.shards, SplitByChromosomeRemote.shards])
            Float splitChromosomeContigLengths = select_first([SplitByChromosome.contig_lengths, SplitByChromosomeRemote.contig_lengths])
            Pair[File, Float] split_chromosomes = (splitChromosomeShards, splitChromosomeContigLengths)
        }
    }

    if (split_into_shards) {
    # if already split into chromosomes, shard further
        if (defined(split_chromosomes)) {
            scatter (chrom_pair in select_first([split_chromosomes])) {
                File chrom_shard = select_first([chrom_pair.left])
                Float chrom_n_records = select_first([chrom_pair.right])
                Int chrom_n_shards = ceil(chrom_n_records / select_first([records_per_shard, 0]))
                call ExecuteScattering as scatterChromosomes {
                    input:
                        vcf_file=chrom_shard,
                        split_vcf_hail_script=split_vcf_hail_script,
                        n_shards=chrom_n_shards,
                        records_per_shard=0,
                        hail_docker=hail_docker,
                        genome_build=genome_build,
                        runtime_attr_override=runtime_attr_split_into_shards
                }
            }
            Array[File] chromosome_shards = flatten(scatterChromosomes.shards)
        }
        
        if (!defined(split_chromosomes)) {
            if (localize_vcf) {
                call ExecuteScattering {
                    input:
                        vcf_file=file,
                        split_vcf_hail_script=split_vcf_hail_script,
                        n_shards=select_first([n_shards]),
                        records_per_shard=select_first([records_per_shard, 0]),
                        hail_docker=hail_docker,
                        genome_build=genome_build,
                        runtime_attr_override=runtime_attr_split_into_shards
                    }
            }
            if (!localize_vcf) {
                String mt_uri = file
                call Helpers.GetHailMTSize as getHailMTSize {
                input:
                    mt_uri=mt_uri,
                    hail_docker=hail_docker
                }
                call ScatterVCFRemote {
                input:
                    vcf_file=mt_uri,
                    input_size=getHailMTSize.mt_size,
                    split_vcf_hail_script=split_vcf_hail_script,
                    n_shards=select_first([n_shards]),
                    records_per_shard=select_first([records_per_shard, 0]),
                    genome_build=genome_build,
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_split_into_shards
                }
            }
        }
    }    

    output {
        Array[File] vcf_shards = select_first([ExecuteScattering.shards, ScatterVCFRemote.shards, chromosome_shards, splitChromosomeShards, [file]])
    }
}   

task GetChromosomeSizes {
    input {
        String vcf_file
        String sv_base_mini_docker
        Boolean has_index
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_disk_gb = 10.0
    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<        
        set -euo pipefail
        if [[ "~{has_index}" == "false" ]]; then
            mkfifo /tmp/token_fifo
            ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
            HTS_AUTH_LOCATION=/tmp/token_fifo tabix --verbosity 3 ~{vcf_file}
        fi;
        export GCS_OAUTH_TOKEN=`/google-cloud-sdk/bin/gcloud auth application-default print-access-token`
        bcftools index -s ~{vcf_file} | cut -f1,3 > contig_lengths.txt
        bcftools query -l ~{vcf_file} | wc -l > n_samples.txt
    >>>

    output {
        Float n_samples = read_lines('n_samples.txt')[0]
        Map[String, Float] contig_lengths = read_map('contig_lengths.txt')
    }
}

task SplitByChromosomeRemote { 
    input {
        String vcf_file
        String chromosome
        String sv_base_mini_docker
        Float input_size
        Boolean has_index
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<        
        set -euo pipefail
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo tabix --verbosity 3 -h ~{vcf_file} ~{chromosome} | bgzip -c > ~{prefix}."~{chromosome}".vcf.gz
        
        tabix ~{prefix}."~{chromosome}".vcf.gz
        # get number of records in chr
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index -n ~{prefix}."~{chromosome}".vcf.gz > contig_length.txt
    >>>

    output {
        File shards = "~{prefix}.~{chromosome}.vcf.gz"
        File shards_idx = "~{prefix}.~{chromosome}.vcf.gz.tbi"
        Float contig_lengths = read_lines('contig_length.txt')[0]
    }
}

task SplitByChromosome { 
    input {
        File vcf_file
        String chromosome
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0
    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        tabix --verbosity 3 ~{vcf_file}
        
        tabix --verbosity 3 -h ~{vcf_file} ~{chromosome} | bgzip -c > ~{prefix}."~{chromosome}".vcf.gz
        
        tabix ~{prefix}."~{chromosome}".vcf.gz
        # get number of records in chr
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index -n ~{prefix}."~{chromosome}".vcf.gz > contig_length.txt
    >>>

    output {
        File shards = "~{prefix}.~{chromosome}.vcf.gz"
        File shards_idx = "~{prefix}.~{chromosome}.vcf.gz.tbi"
        Float contig_lengths = read_lines('contig_length.txt')[0]
    }
}

task ExecuteScattering {
    input {
        File vcf_file
        Int n_shards
        Int records_per_shard
        String split_vcf_hail_script
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB") 
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        curl  ~{split_vcf_hail_script} > split_vcf.py
        python3 split_vcf.py ~{vcf_file} ~{n_shards} ~{records_per_shard} ~{prefix} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} ~{genome_build}
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }
}

task ScatterVCFRemote {
    input {
        String vcf_file
        Float input_size
        Int n_shards
        Int records_per_shard
        String split_vcf_hail_script
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    String prefix = if sub(vcf_file, '.mt', '')!=vcf_file then basename(vcf_file, '.mt') else basename(vcf_file, ".vcf.gz")

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        curl  ~{split_vcf_hail_script} > split_vcf.py
        python3 split_vcf.py ~{vcf_file} ~{n_shards} ~{records_per_shard} ~{prefix} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} ~{genome_build}
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>
    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }
}
