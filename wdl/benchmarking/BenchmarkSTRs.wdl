version 1.0

import "../utils/Structs.wdl"

workflow BenchmarkSTRs {
    input {
        Array[String] sample_ids
        Array[File] trgt_vcfs
        Array[File] trgt_vcf_indices
        File vamos_vcf
        File vamos_vcf_index
        File ref_fa
        File ref_fai
        Array[String] contigs
        
        String output_prefix
        Boolean include_all_regions = false
        
        String benchmark_strs_docker
        
        RuntimeAttr? runtime_attr_subset_vamos
        RuntimeAttr? runtime_attr_run_benchmark
        RuntimeAttr? runtime_attr_aggregate_and_plot
    }
    
    scatter (idx in range(length(sample_ids))) {
        String sample_id = sample_ids[idx]
        File trgt_vcf = trgt_vcfs[idx]
        File trgt_vcf_index = trgt_vcf_indices[idx]
        
        call SubsetAndReheaderVamos {
            input:
                vamos_vcf = vamos_vcf,
                vamos_vcf_index = vamos_vcf_index,
                sample_id = sample_id,
                docker = benchmark_strs_docker,
                runtime_attr_override = runtime_attr_subset_vamos
        }
        
        call RunBenchmark {
            input:
                vamos_vcf = SubsetAndReheaderVamos.subsetted_vcf,
                vamos_vcf_index = SubsetAndReheaderVamos.subsetted_vcf_index,
                trgt_vcf = trgt_vcf,
                trgt_vcf_index = trgt_vcf_index,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                sample_id = sample_id,
                contigs = contigs,
                docker = benchmark_strs_docker,
                runtime_attr_override = runtime_attr_run_benchmark
        }
    }
    
    call AggregateAndPlot {
        input:
            stats_tsvs = RunBenchmark.stats_tsv,
            sample_ids = sample_ids,
            output_prefix = output_prefix,
            include_all_regions = include_all_regions,
            docker = benchmark_strs_docker,
            runtime_attr_override = runtime_attr_aggregate_and_plot
    }
    
    output {
        Array[File] per_sample_stats = RunBenchmark.stats_tsv
        File aggregated_match_data_non_ref = AggregateAndPlot.match_data_non_ref
        File genotype_concordance_matrix_non_ref = AggregateAndPlot.genotype_concordance_matrix_non_ref
        File similarity_plot_non_ref = AggregateAndPlot.similarity_plot_non_ref
        File edit_distance_plot_non_ref = AggregateAndPlot.edit_distance_plot_non_ref
        File length_difference_plot_non_ref = AggregateAndPlot.length_difference_plot_non_ref
        File length_diff_vs_locus_size_non_ref = AggregateAndPlot.length_diff_vs_locus_size_non_ref
        
        File? aggregated_match_data_all = AggregateAndPlot.match_data_all
        File? genotype_concordance_matrix_all = AggregateAndPlot.genotype_concordance_matrix_all
        File? similarity_plot_all = AggregateAndPlot.similarity_plot_all
        File? edit_distance_plot_all = AggregateAndPlot.edit_distance_plot_all
        File? length_difference_plot_all = AggregateAndPlot.length_difference_plot_all
        File? length_diff_vs_locus_size_all = AggregateAndPlot.length_diff_vs_locus_size_all
        File? edit_distance_to_reference_all = AggregateAndPlot.edit_distance_to_reference_all
        File? length_difference_to_reference_all = AggregateAndPlot.length_difference_to_reference_all
    }
}

task SubsetAndReheaderVamos {
    input {
        File vamos_vcf
        File vamos_vcf_index
        String sample_id
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size([vamos_vcf, vamos_vcf_index], "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 3.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 1,
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }
    
    String output_vcf = "~{sample_id}.vamos.subsetted.vcf.gz"
    
    command <<<
        set -euo pipefail
        
        # Get list of samples in VCF
        bcftools query -l ~{vamos_vcf} > samples.txt
        
        # Find the sample matching pattern: _SAMPLE_ID_
        # The sample ID should be embedded within underscores
        matched_sample=$(grep -E "_~{sample_id}_" samples.txt | head -n1)
        
        if [ -z "$matched_sample" ]; then
            echo "ERROR: Could not find sample matching pattern _~{sample_id}_ in VCF"
            exit 1
        fi
        
        echo "Found matching sample: $matched_sample"
        
        # Subset to the matched sample and GRCh38_to_GRCh38/GRCh38_to_GRCh38
        bcftools view -s "$matched_sample,GRCh38_to_GRCh38/GRCh38_to_GRCh38" \
            ~{vamos_vcf} \
            -O z \
            -o temp.vcf.gz
        
        bcftools index -t temp.vcf.gz
        
        # Create reheader file to rename the sample
        echo "$matched_sample ~{sample_id}" > reheader.txt
        echo "GRCh38_to_GRCh38/GRCh38_to_GRCh38 GRCh38_to_GRCh38/GRCh38_to_GRCh38" >> reheader.txt
        
        # Reheader the VCF
        bcftools reheader -s reheader.txt temp.vcf.gz -o ~{output_vcf}
        bcftools index -t ~{output_vcf}
    >>>
    
    output {
        File subsetted_vcf = output_vcf
        File subsetted_vcf_index = "~{output_vcf}.tbi"
    }
}

task RunBenchmark {
    input {
        File vamos_vcf
        File vamos_vcf_index
        File trgt_vcf
        File trgt_vcf_index
        File ref_fa
        File ref_fai
        String sample_id
        Array[String] contigs
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size([vamos_vcf, vamos_vcf_index, trgt_vcf, trgt_vcf_index, ref_fa, ref_fai], "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 2,
        preemptible_tries: 1,
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        
        python3 /opt/gnomad-lr/scripts/benchmark/benchmark_str_concordance.py \
            --vamos-vcf ~{vamos_vcf} \
            --trgt-vcf ~{trgt_vcf} \
            --ref-fa ~{ref_fa} \
            --contigs ~{sep="," contigs} \
            --output-stats ~{sample_id}.stats.tsv \
            --sample-id ~{sample_id}
    >>>
    
    output {
        File stats_tsv = "~{sample_id}.stats.tsv"
    }
}

task AggregateAndPlot {
    input {
        Array[File] stats_tsvs
        Array[String] sample_ids
        String output_prefix
        Boolean include_all_regions
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size(stats_tsvs, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 3.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 2,
        preemptible_tries: 1,
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }
    
    command <<<
        set -euo pipefail
        
        # Create file with all stats TSVs
        cat ~{write_lines(stats_tsvs)} > stats_files.txt
        
        python3 /opt/gnomad-lr/scripts/benchmark/aggregate_str_benchmarks.py \
            --stats-files stats_files.txt \
            --output-prefix ~{output_prefix} \
            --include-all-regions ~{if include_all_regions then "true" else "false"}
    >>>
    
    output {
        File match_data_non_ref = "~{output_prefix}.non_ref_regions.match_data.tsv"
        File genotype_concordance_matrix_non_ref = "~{output_prefix}.non_ref_regions.genotype_concordance_matrix.tsv"
        File similarity_plot_non_ref = "~{output_prefix}.non_ref_regions.similarity_across_callers.png"
        File edit_distance_plot_non_ref = "~{output_prefix}.non_ref_regions.edit_distance_across_callers.png"
        File length_difference_plot_non_ref = "~{output_prefix}.non_ref_regions.length_difference_across_callers.png"
        File length_diff_vs_locus_size_non_ref = "~{output_prefix}.non_ref_regions.length_difference_vs_locus_size.png"
        
        File? match_data_all = if include_all_regions then "~{output_prefix}.all_regions.match_data.tsv" else ""
        File? genotype_concordance_matrix_all = if include_all_regions then "~{output_prefix}.all_regions.genotype_concordance_matrix.tsv" else ""
        File? similarity_plot_all = if include_all_regions then "~{output_prefix}.all_regions.similarity_across_callers.png" else ""
        File? edit_distance_plot_all = if include_all_regions then "~{output_prefix}.all_regions.edit_distance_across_callers.png" else ""
        File? length_difference_plot_all = if include_all_regions then "~{output_prefix}.all_regions.length_difference_across_callers.png" else ""
        File? length_diff_vs_locus_size_all = if include_all_regions then "~{output_prefix}.all_regions.length_difference_vs_locus_size.png" else ""
        File? edit_distance_to_reference_all = if include_all_regions then "~{output_prefix}.all_regions.edit_distance_to_reference.png" else ""
        File? length_difference_to_reference_all = if include_all_regions then "~{output_prefix}.all_regions.length_difference_to_reference.png" else ""
    }
}
