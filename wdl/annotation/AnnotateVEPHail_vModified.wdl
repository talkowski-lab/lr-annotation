version 1.0
    
import "../utils/Structs.wdl"
import "../utils/ScatterVCF.wdl" as ScatterVCF
import "../utils/Helpers.wdl" as Helpers

workflow AnnotateVEPHail_vModified {
    input {
        File vcf
        File vcf_idx

        String cohort_prefix
        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/vep/split_vcf_hail.py"
        String vep_annotate_hail_python_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/vep/vep_annotate_hail_vModified.py"
        String genome_build = "GRCh38"
        Boolean split_by_chromosome
        Boolean split_into_shards

        File top_level_fa
        File ref_vep_cache

        String annotate_vep_hail_docker
        String hail_docker
        String sv_base_mini_docker
        
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
                docker=annotate_vep_hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_vep_annotate
        }
    }
    
    call Helpers.ConcatTsvs {
        input:
            tsvs=VepAnnotate.vep_tsv_file,
            prefix=cohort_prefix + ".vep_annotations",
            docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_combine_vcfs
    }

    output {
        File annotations_tsv_vep = ConcatTsvs.concatenated_tsv
    }
}   

task VepAnnotate {
    input {
        File vcf
        File top_level_fa
        File ref_vep_cache
        String genome_build
        String vep_annotate_hail_python_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB") + size(ref_vep_cache, "GB")) + 20,
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

    String filename = basename(vcf)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf, ".vcf.gz") else basename(vcf, ".vcf.bgz")

    command <<<
        set -euo pipefail
        
        mkdir -p cache
        tar xzf ~{ref_vep_cache} -C ./cache

        dir_cache_path=$(pwd)/cache

        echo '{"command": [
        "/opt/vep/src/ensembl-vep/vep",
            "--format", "vcf",
            "__OUTPUT_FORMAT_FLAG__",
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
            "-o", "STDOUT"
        ],
        "env": {},
        "vep_json_schema": "Struct{allele_string:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,id:String,phenotype_or_disease:Int32,pubmed:Array[Int32],somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,flags:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,mane_select:String,mane_plus_clinical:String,minimised:Int32,mirna:Array[String],polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,uniprot_isoform:Array[String],variant_allele:String}],variant_class:String}"
        }' > vep_config.json

        curl ~{vep_annotate_hail_python_script} > vep_annotate.py

        python3 vep_annotate.py \
            -i ~{vcf} \
            -o ~{prefix}.vep.tsv \
            --cores ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} \
            --build ~{genome_build}
                
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File vep_tsv_file = "~{prefix}.vep.tsv"
        File hail_log = "hail_log.txt"
    }
}
