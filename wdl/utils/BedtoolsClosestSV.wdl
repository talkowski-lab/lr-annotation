version 1.0

import "Helpers.wdl"
import "Structs.wdl"

workflow BedtoolsClosestSV {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx
        String prefix

        Int min_sv_length_eval
        Int min_sv_length_truth
        String type_field_eval = "SVTYPE"
        String length_field_eval = "SVLEN"
        
        String sv_pipeline_docker
        String utils_docker
        
        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_match
    }

    call Helpers.SubsetVcfByLength as SubsetEval {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval_idx,
            length_field = length_field_eval,
            min_length = min_sv_length_eval,
            prefix = "~{prefix}.subset_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_eval
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = SubsetEval.subset_vcf,
            vcf_idx = SubsetEval.subset_vcf_idx,
            type_field = type_field_eval,
            length_field = length_field_eval,
            prefix = "~{prefix}.eval.symbolic",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call Helpers.SubsetVcfByLength as SubsetTruth {
        input:
            vcf = vcf_sv_truth,
            vcf_idx = vcf_sv_truth_idx,
            length_field = "SVLEN",
            min_length = min_sv_length_truth,
            prefix = "~{prefix}.subset_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_truth
    }

    call AnnotateExternalAFMatches {
        input:
            vcf_eval = ConvertToSymbolic.processed_vcf,
            vcf_eval_idx = ConvertToSymbolic.processed_vcf_idx,
            vcf_truth = SubsetTruth.subset_vcf,
            vcf_truth_idx = SubsetTruth.subset_vcf_idx,
            prefix = prefix,
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_match
    }

    output {
        File closest_bed = AnnotateExternalAFMatches.closest_bed
        File annotation_tsv = AnnotateExternalAFMatches.annotation_tsv
    }
}

task AnnotateExternalAFMatches {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export LC_ALL=C

        split_vcf_to_beds() {
            local vcf="$1"
            local prefix="$2"
            local include_af="$3"

            if [[ "$include_af" == "true" ]]; then
                svtk vcf2bed -i SVTYPE -i SVLEN -i AF "$vcf" tmp.${prefix}.bed
                awk 'BEGIN{FS=OFS="\t"} NR == 1 {print $1,$2,$3,$4,"svtype",$7,$8,$9; next} {af=$9; if (af == "") af="."; print $1,$2,$3,$4,$7,$7,$8,af}' tmp.${prefix}.bed > ${prefix}.bed
            else
                svtk vcf2bed -i SVTYPE -i SVLEN "$vcf" tmp.${prefix}.bed
                cut -f1-4,7-8 tmp.${prefix}.bed > ${prefix}.bed
            fi

            head -1 ${prefix}.bed > ${prefix}.header
            type_col=5
            if [[ "$include_af" == "true" ]]; then
                type_col=6
            fi
            cat ${prefix}.header <(awk -v c="$type_col" 'toupper($c) == "DEL"' ${prefix}.bed) > ${prefix}.DEL.bed
            cat ${prefix}.header <(awk -v c="$type_col" 'toupper($c) == "DUP"' ${prefix}.bed) > ${prefix}.DUP.bed
            cat ${prefix}.header <(awk -v c="$type_col" 'toupper($c) ~ /^(INS|INS:ME|INS:ME:ALU|INS:ME:LINE1|INS:ME:SVA|ALU|LINE1|SVA|HERVK)$/' ${prefix}.bed) > ${prefix}.INS.bed
            cat ${prefix}.header <(awk -v c="$type_col" 'toupper($c) == "INV" || toupper($c) == "CPX"' ${prefix}.bed) > ${prefix}.INV.bed
            cat ${prefix}.header <(awk -v c="$type_col" 'toupper($c) == "BND" || toupper($c) == "CTX"' ${prefix}.bed) > ${prefix}.BND.bed
        }

        run_closest() {
            local svtype="$1"
            local selector="$2"
            local query_bed="eval.${svtype}.bed"
            local truth_bed="truth.${svtype}.bed"
            local closest_bed="${svtype}.closest.bed"
            local comparison="${svtype}.comparison"

            paste <(head -1 "$query_bed") <(head -1 "$truth_bed") \
                | sed -e "s/#//g" \
                | awk 'BEGIN{OFS="\t"} {print $0,"overlap"}' \
                > "$closest_bed"
            expected_columns=$(awk -F'\t' 'NR == 1 {print NF}' "$closest_bed")

            tail -n +2 "$query_bed" | awk 'NF > 0' | sort -k1,1 -k2,2n > query.body.bed
            tail -n +2 "$truth_bed" | awk 'NF > 0' | sort -k1,1 -k2,2n > truth.body.bed

            if [[ -s query.body.bed && -s truth.body.bed ]]; then
                bedtools closest \
                    -wo \
                    -a query.body.bed \
                    -b truth.body.bed \
                | awk -v expected_columns="$expected_columns" 'BEGIN{FS=OFS="\t"} NF == expected_columns {print} NF != expected_columns {print "Skipping malformed bedtools closest row with " NF " columns: " $0 > "/dev/stderr"}' \
                >> "$closest_bed"
            fi

            Rscript "$selector" \
                -i "$closest_bed" \
                -o "$comparison" \
                -p pop
        }

        split_vcf_to_beds ~{vcf_eval} eval false
        split_vcf_to_beds ~{vcf_truth} truth true
        echo "ALL" > pop

        run_closest DEL /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R
        run_closest DUP /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R
        run_closest INS /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R
        run_closest INV /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R
        run_closest BND /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R

        awk 'FNR == 1 && NR != 1 {next} {print}' \
            DEL.comparison \
            DUP.comparison \
            INS.comparison \
            INV.comparison \
            BND.comparison \
            > ~{prefix}.comparison.tsv

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{vcf_eval} \
        | sort -k5,5 > eval_info.tsv

        awk -F'\t' 'BEGIN{OFS="\t"} $1 != "query_svid" && $1 != "" && $2 != "" {
            truth_af = ($3 == "" ? "." : $3)
            print $1, $2, truth_af
        }' ~{prefix}.comparison.tsv \
            | sort -k1,1 > matched_ids.tsv

        join \
            -1 5 \
            -2 1 \
            -t $'\t' \
            eval_info.tsv \
            matched_ids.tsv \
            > joined.tsv

        awk -F'\t' 'BEGIN{OFS="\t"} {
            print $2, $3, $4, $5, $1, "BEDTOOLS_CLOSEST", $6, "SV", $7
        }' joined.tsv > ~{prefix}.bedtools_matched.tsv
    >>>

    output {
        File closest_bed = "~{prefix}.comparison.tsv"
        File annotation_tsv = "~{prefix}.bedtools_matched.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
