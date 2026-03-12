version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateGQMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Array[String] gq_fields
        Array[Array[Int]] gq_bins
        Array[String] gq_variant_filters
        Array[Boolean] gq_larger_field

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call AnnotateGQMetricsTask {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.gq_annotated",
                gq_fields = gq_fields,
                gq_bins = gq_bins,
                gq_variant_filters = gq_variant_filters,
                gq_larger_field = gq_larger_field,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateGQMetricsTask.annotated_vcf,
            vcf_idxs = AnnotateGQMetricsTask.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.gq_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File gq_annotated_vcf = ConcatVcfs.concat_vcf
        File gq_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateGQMetricsTask {
    input {
        File vcf
        File vcf_idx
        String prefix

        Array[String] gq_fields
        Array[Array[Int]] gq_bins
        Array[String] gq_variant_filters
        Array[Boolean] gq_larger_field

        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # 1. Provide array variables to bash via file writes
        cat ~{write_lines(gq_fields)} > fields.txt
        cat ~{write_lines(gq_variant_filters)} > filters.txt
        cat ~{write_lines(gq_larger_field)} > larger.txt
        cat ~{write_json(gq_bins)} > bins.json

        # 2. Setup standard boilerplate python script to calculate metric aggregates per field
        cat << 'EOF' > process_field.py
        import sys
        import json

        field = sys.argv[1]
        larger_flag = sys.argv[2].lower() == 'true'
        bin_idx = int(sys.argv[3])

        with open("bins.json") as f:
            all_bins = json.load(f)
        bins = all_bins[bin_idx]

        # Define VCF-compliant Headers (Number=A means it expects a comma-separated list of values, one per ALT)
        header = f'##INFO=<ID={field}_hist_all_bin_freq,Number=A,Type=String,Description="Histogram for {field} all individuals; bin edges are: 0|{"|".join(map(str, bins))}">\n'
        header += f'##INFO=<ID={field}_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for {field} in heterozygous/alt individuals; bin edges are: 0|{"|".join(map(str, bins))}">\n'
        if larger_flag:
            header += f'##INFO=<ID={field}_hist_all_n_larger,Number=A,Type=Integer,Description="Count of all genotypes with {field} larger than {bins[-1]}">\n'
            header += f'##INFO=<ID={field}_hist_alt_n_larger,Number=A,Type=Integer,Description="Count of alt genotypes with {field} larger than {bins[-1]}">\n'

        with open(f"header_{bin_idx}.txt", "w") as f:
            f.write(header)

        with open(f"cols_{bin_idx}.txt", "w") as f:
            cols = f"CHROM,POS,REF,ALT,INFO/{field}_hist_all_bin_freq,INFO/{field}_hist_alt_bin_freq"
            if larger_flag:
                 cols += f",INFO/{field}_hist_all_n_larger,INFO/{field}_hist_alt_n_larger"
            f.write(cols)

        out_file = open(f"annots_{bin_idx}.tsv", "w")
        for line in sys.stdin:
            parts = line.strip('\n').split('\t')
            if len(parts) < 4: continue
            chrom, pos, ref, alt = parts[:4]
            samples = parts[4:]
            
            counts_all = [0] * (len(bins) + 1)
            counts_alt = [0] * (len(bins) + 1)
            
            for sample in samples:
                if '|' not in sample: continue
                gt_str, val_str = sample.split('|', 1)
                
                # Check missing/uncalled
                if val_str == '.' or val_str == '': continue
                try:
                    val = float(val_str)
                except ValueError:
                    continue
                
                # Exclude purely REF or purely missing GTs
                is_alt = False
                for a in gt_str.replace('|', '/').split('/'):
                    if a != '.' and a != '0':
                        is_alt = True
                        break
                
                # Determine bin
                b_idx = len(bins)
                for j, b in enumerate(bins):
                    if val <= b:
                        b_idx = j
                        break
                
                counts_all[b_idx] += 1
                if is_alt:
                    counts_alt[b_idx] += 1
                    
            all_str = "|".join(map(str, counts_all))
            alt_str = "|".join(map(str, counts_alt))
            
            # Replicate values per ALT allele sequentially to appease the Number=A configuration requested.
            num_alts = len(alt.split(','))
            out_all = ",".join([all_str] * num_alts)
            out_alt = ",".join([alt_str] * num_alts)
            
            out_row = [chrom, pos, ref, alt, out_all, out_alt]
            
            # Evaluate extra requirements dynamically
            if larger_flag:
                val_all_larger = str(counts_all[-1])
                val_alt_larger = str(counts_alt[-1])
                out_row.extend([",".join([val_all_larger] * num_alts), ",".join([val_alt_larger] * num_alts)])
                
            out_file.write("\t".join(out_row) + "\n")
        out_file.close()
        EOF

        cp ~{vcf} current.vcf.gz
        cp ~{vcf_idx} current.vcf.gz.tbi

        NUM_FIELDS=$(wc -l < fields.txt)

        # 3. Iterate through arrays mapping Python evaluation sequentially
        for i in $(seq 0 $((NUM_FIELDS - 1))); do
            FIELD=$(sed -n "$((i+1))p" fields.txt)
            FILTER=$(sed -n "$((i+1))p" filters.txt)
            LARGER=$(sed -n "$((i+1))p" larger.txt)
            
            FILTER_ARGS=""
            if [ -n "$FILTER" ] && [ "$FILTER" != "." ] && [ "$FILTER" != "None" ]; then
                FILTER_ARGS="-i '${FILTER}'"
            fi
            
            # Conditionally supply bcftools expression filtering
            if [ -n "$FILTER_ARGS" ]; then
                bcftools query $FILTER_ARGS -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT|%${FIELD}]\n" current.vcf.gz | \
                    python process_field.py "$FIELD" "$LARGER" "$i"
            else
                bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT|%${FIELD}]\n" current.vcf.gz | \
                    python process_field.py "$FIELD" "$LARGER" "$i"
            fi
                
            bgzip -c annots_${i}.tsv > annots_${i}.tsv.gz
            tabix -s1 -b2 -e2 annots_${i}.tsv.gz
            
            COLS=$(cat cols_${i}.txt)
            
            # Incrementally annotate file 
            bcftools annotate -a annots_${i}.tsv.gz -h header_${i}.txt -c "$COLS" \
                -Oz -o next.vcf.gz current.vcf.gz
            
            tabix -p vcf next.vcf.gz
            mv next.vcf.gz current.vcf.gz
            mv next.vcf.gz.tbi current.vcf.gz.tbi
        done

        mv current.vcf.gz ~{prefix}.vcf.gz
        mv current.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
