version 1.0
    
import "general/Structs.wdl"
import "general/ScatterVCF.wdl" as ScatterVCF
import "MergeSplitVCF.wdl" as MergeSplitVCF

workflow AnnotateVEPHail {
    input {
        File vcf
        File vcf_index

        File top_level_fa
        File ref_vep_cache

        String cohort_prefix
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        
        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/annotate_vep_hail/split_vcf_hail.py"
        String vep_annotate_hail_python_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/annotate_vep_hail/vep_annotate_hail.py"

        String genome_build="GRCh38"
        Boolean split_by_chromosome
        Boolean split_into_shards 

        RuntimeAttr? runtime_attr_split_by_chr
        RuntimeAttr? runtime_attr_split_into_shards
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_combine_vcfs
    }

    call ScatterVCF.ScatterVCF {
        input:
            file=vcf,
            has_index=true,
            split_vcf_hail_script=split_vcf_hail_script,
            cohort_prefix=cohort_prefix,
            genome_build=genome_build,
            split_by_chromosome=split_by_chromosome,
            split_into_shards=split_into_shards,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_split_by_chr=runtime_attr_split_by_chr,
            runtime_attr_split_into_shards=runtime_attr_split_into_shards
    }

    scatter (vcf_shard in ScatterVCF.vcf_shards) {
        call VepAnnotate {
            input:
                vcf=vcf_shard,
                vep_annotate_hail_python_script=vep_annotate_hail_python_script,
                top_level_fa=top_level_fa,
                ref_vep_cache=ref_vep_cache,
                vep_hail_docker=vep_hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_vep_annotate
        }
    }
    
    call MergeSplitVCF.CombineVCFs {
        input:
            vcf_files=VepAnnotate.vep_vcf_file,
            vcf_indices=VepAnnotate.vep_vcf_idx,
            naive=true,
            allow_overlaps=false,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + ".final",
            sort_after_merge=true,
            runtime_attr_override=runtime_attr_combine_vcfs
    }

    output {
        File vep_annotated_vcf = CombineVCFs.combined_vcf
        File vep_annotated_vcf_index = CombineVCFs.combined_vcf_index
    }
}   

task VepAnnotate {
    input {
        File vcf
        File top_level_fa
        File ref_vep_cache
        String genome_build
        String vep_annotate_hail_python_script
        String vep_hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB") + size(ref_vep_cache, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr default_attr = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    String filename = basename(vcf)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf, ".vcf.gz") else basename(vcf, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.vep.vcf.bgz"

    command <<<
        set -euo pipefail
        
        mkdir -p cache
        tar xzf ~{ref_vep_cache} -C ./cache

        dir_cache_path=$(pwd)/cache

        echo '{"command": [
        "/opt/vep/src/ensembl-vep/vep",
            "--format", "vcf",
            "--json",
            "--everything",
            "--allele_number",
            "--no_stats",
            "--cache", 
            "--offline",
            "--minimal",
            "--assembly", "~{genome_build}",
            "--merged",
            "--fasta", "~{top_level_fa}",
            "--dir_cache", "'$dir_cache_path'",
            "--vcf_info_field", "vep",
            "-o", "STDOUT"
        ],
        "env": {},
        "vep_json_schema": "Struct{allele_string:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,id:String,phenotype_or_disease:Int32,pubmed:Array[Int32],somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,flags:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,mane_select:String,mane_plus_clinical:String,minimised:Int32,mirna:Array[String],polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,uniprot_isoform:Array[String],variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        curl ~{vep_annotate_hail_python_script} > vep_annotate.py

        proj_id=$(gcloud config get-value project)

        python3 vep_annotate.py \
            -i ~{vcf} \
            -o ~{vep_annotated_vcf_name} \
            --cores ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} \
            --mem ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} \
            --build ~{genome_build} \
            --project-id $proj_id
        
        cp $(ls . | grep hail*.log) hail_log.txt

        bcftools index -t ~{vep_annotated_vcf_name}
    >>>

    output {
        File vep_vcf_file = vep_annotated_vcf_name
        File vep_vcf_idx = vep_annotated_vcf_name+".tbi"
        File hail_log = "hail_log.txt"
    }
}
