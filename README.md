# Long-Read Annotation
This repository serves as a home for all scripts, workflows and processes for annotating long-read callsets.



## Cohort
- HPRC.
	- 232 total samples.
	- [Metadata File](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
		- 234 total samples.
		- Additionally includes GRC38 and CHM13.
		- Fabio's table and the DeepVariant callset also include HG03492.
- HGSVC.
	- 65 total samples.
	- [Metadata File](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
		- 67 total samples.
		- Renames NA21487 to GM21487.
		- Duplicates NA19129 (also includes GM19129) and NA20355 (also includes GM20355).
		- Additionally includes GM19320.
		- Misses NA24385 (HG002 in Terra).
- HPRC & HGSVC Overlapping Samples: HG002, HG00733, HG02818, NA19036, NA19240.
- All of Us Phase 1.
	- 1027 total samples.
	- Metadata file extracted from VCF - all samples are unrelated and of African ancestry.


## References
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.
- `dbsnp_vcf` (per-contig in data table): [Build 156 variants for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GCF_000001405.40.vcf.gz) from [dbSNP](https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz).
- `dbsnp_vcf_idx` (per-contig in data table): Index for `dbsnp_vcf`.
- `dbvar_vcf` (per-contig in data table): [Variants for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/GRCh38.variant_call.all.vcf.gz) from [dbVaR](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz).
- `dbvar_vcf_idx` (per-contig in data table): Index for `dbvar_vcf`.
- `exons_bed`: [Loci for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `fix_variant_collisions_java`: [Script for phasing](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/FixVariantCollisions-9daf9b2.java) from the DSP Long-Read SV team.
- `genetic_maps_tsv`: [Maps for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/genetic_map_b38.tsv) from the [GLIMPSE references](https://github.com/odelaneau/GLIMPSE/tree/master/maps/genetic_maps.b38).
- `gnomad_sv_constraint`: Genic constraint results for `gnomad_sv_vcf`.
- `gnomad_sv_vcf`: [gnomAD V4 structural variants](gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz) from the [gnomAD downloads](https://gnomad.broadinstitute.org/data#v4). 
- `gnomad_sv_vcf_idx`: Index for `gnomad_sv_vcf`.
- `gnomad_tr_json`: [gnomAD V4 tandem repeats](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variant_catalog_without_offtargets.GRCh38.json) from the [str-analysis references](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json).
- `gnomad_vcf` (per-contig in data table): [gnomAD V4 short variants](gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/) from the [gnomAD downloads](https://gnomad.broadinstitute.org/data#v4). 
- `gnomad_vcf_idx` (per-contig in data table): Indexes for `gnomad_vcfs`.
- `insilico_cadd_ht`: [CADD scores for GRCh38](gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.cadd.grch38.ht) from gnomAD V4.
- `insilico_pangolin_ht`: [Pangolin scores for GRCh38](gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.pangolin.grch38.ht) from gnomAD V4.
- `insilico_phylop_ht`: [PhyloP scores for GRCh38](gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.phylop.grch38.ht) from gnomAD V4.
- `insilico_revel_ht`: [REVEL scores for GRCh38](gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.revel.grch38.ht) from gnomAD V4.
- `insilico_spliceai_ht`: [SpliceAI scores for GRCh38](gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.spliceai.grch38.ht) from gnomAD V4.
- `mei_catalog`: [Loci for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/meis_hg38.tsv) is a combination of ALU and LINE loci from [RepeatMasker](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema) and SVA loci from [van Bree et al](https://genome.cshlp.org/content/suppl/2022/03/24/gr.275515.121.DC1).
- `mei_fa`: [Sequences for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).
- `mei_fa_indices`: BWA and Minimap indices for `mei_fa`.
- `noncoding_bed`: [Panel for GRCh38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `par_bed`: [Panel for GRCh38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `ploidy_bed_female`: [Panel for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/grch38_female.bed) from the [Kanpig references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `ploidy_bed_male`: [Panel for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/grch38_male.bed) from the [Kanpig references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `ref_dict`: Sequence dictionary for `ref_fa`.
- `ref_fa`: [Sequences for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/hg38.no_alt.fa) from the [PAV references](https://github.com/ACEnglish/kanpig/tree/main/ploidy_beds).
- `ref_fai`: Index for `ref_fa`.
- `ref_fa_bwa_indices`: BWA indices for `ref_fa`.
- `ref_fa_indices`: BWA and Minimap indices for `ref_fa`.
- `ref_vep_cache`: [Cache for v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from the [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/), with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).
- `repeat_masker_bed`: [Loci for GRCh38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/GenomicContext/hg38.RM.sorted.merged.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `repeats_bed`: [Loci for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).
- `repeat_catalog_trgt`: [Catalog for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [the TR Catalog references](https://github.com/broadinstitute/tandem-repeat-catalog/releases), as used in All of Us Phase.
- `seg_dup_bed`: [Loci for GRCh38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/GenomicContext/hg38.SegDup.sorted.merged.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `simple_repeat_bed`: [Loci for GRCh38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/GenomicContext/hg38.SimpRep.sorted.merged.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `seqrepo_tar`: [Sequence repository v2024-12-20 for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/seqrepo_2024-12-20.tar.gz) from the [latest release of seqrepo](https://github.com/biocommons/biocommons.seqrepo).
- `truth_hgsvc_ped`: Sample pedigree file from the [HGSVC FTP server](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/).
- `truth_hgsvc_vcf`: [HGSVC3 release](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/), combined across the specific VCFs for SNVs, indels and SVs, from the [HGSVC FTP server](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/).
- `truth_hgsvc_vcf_idx`: Index for `truth_hgsvc_vcf`.
- `truth_hprc_ped`: Sample pedigree file from the [HPRC AWS bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=).
- `truth_hprc_vcf`: [HPRC2 release](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.wave.vcf.gz) from the [HPRC AWS bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=).
- `truth_hprc_vcf_idx`: Index for `truth_hprc_vcf`.
- `vntr_bed`: [Loci for GRCh38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).



## Annotation Workflows
### [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow leverages [AnnotateVcf](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/AnnotateVcf.wdl) from the GATK-SV pipeline in order to annotate internal allele frequencies based on sample sexes and ancestries. It runs on all variants in the input VCF, including SVs.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `File sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `File ped_file`: Six column file containing the cohort pedigree, with specifications described in [this article](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
- `File lps_tsv`: LPS TSV file output by the [`trgt-lps` tool](https://github.com/PacificBiosciences/trgt-lps), which contains information regarding the longest polymer sequence within each sample for each loci genotyped by TRGT.
- `Array[String]? strip_info_fields`: Comma-separated list of `INFO` fields to remove from the input VCF prior to annotation.
- `Int records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File par_bed`: From [references](#references).

Outputs:
- `af_annotated_vcf`: Annotated VCF.
- `af_annotated_vcf_idx`: Index for annotated VCF.


### [AnnotateAgeMetrics](wdl/annotation/AnnotateAgeMetrics.wdl)
This workflow computes the age distribution of carriers for every variant in the input VCF. For each sample it derives an age from a date-of-birth table relative to a fixed reference date, then tabulates the number of heterozygous and homozygous carriers of each allele that fall into a set of user-defined age bins, along with overflow `smaller` and `larger` bins for ages outside the configured range. It emits a TSV of these per-allele age-bin counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `File age_data`: CSV file with `person_id` and `date_of_birth` columns, used to derive each sample's age.
- `Array[Int] age_bins`: Age-bin edges, in years, into which carrier ages are binned.
- `String reference_date`: Reference date (`YYYY-MM-DD`) against which each sample's age is computed.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_age`: TSV of per-allele carrier counts across the age bins.


### [AnnotateCallsetOverlap](wdl/annotation/AnnotateCallsetOverlap.wdl)
This workflow ingests an evaluation VCF and two truth VCFs — one of SNVs & indels and one of SVs — and finds matching variants across them in order to compare the AF & VEP annotations of the matched pairs. This serves as a degree of benchmarking, as it ensures that annotations applied to a larger cohort (e.g. gnomAD) are in line with those we annotate. It also enables the identification of variants that are outliers relative to existing cohorts by pulling out those with a large amount of discordance in their annotation across the callsets.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs:
1. Exact match across CHROM, POS, REF and ALT.
2. Truvari match with overlap percentages of 90%, 70% and 50%.
3. Matching based on `bedtools closest`, finetuned for SVs. Here the evaluation and truth variants are split by type and converted to a symbolic representation, after which separate `bedtools closest` passes are run — one tuned for deletions and duplications via reciprocal positional overlap, and one tuned for insertions via breakpoint proximity — so that each evaluation variant is paired with the nearest same-type truth variant above the per-callset minimum SV-length thresholds.

> **Note:** This workflow considers the DUP region for INS annotated as DUP.

Inputs:
- `File vcf_eval`: VCF whose annotations are being evaluated.
- `File vcf_eval_idx`: Index for `vcf_eval`.
- `File vcf_truth`: VCF containing SNVs & indels to evaluate against.
- `File vcf_truth_idx`: Index for `vcf_truth`.
- `File vcf_sv_truth`: VCF containing SVs to evaluate against.
- `File vcf_sv_truth_idx`: Index for `vcf_sv_truth`.
- `Array[String] contigs`: Contigs to evaluate.
- `Int min_sv_length_eval_truvari`: Minimum length for an evaluation variant to enter the Truvari matching round.
- `Int min_sv_length_truth_truvari`: Minimum length for a truth variant to enter the Truvari matching round.
- `Int min_sv_length_eval_bedtools_closest`: Minimum length for an evaluation variant to enter the `bedtools closest` matching round.
- `Int min_sv_length_truth_bedtools_closest`: Minimum length for a truth variant to enter the `bedtools closest` matching round.
- `Boolean compare_annotations`: Whether to compute the downstream AF & VEP annotation comparison summaries (default `true`).
- `Boolean do_exact`: Whether to run the exact matching round (default `true`).
- `Boolean do_truvari`: Whether to run the Truvari matching round (default `true`).
- `Boolean do_bedtools_closest`: Whether to run the `bedtools closest` matching round (default `true`).
- `String type_field_eval`: INFO field in the evaluation VCF giving each variant's allele type (default `allele_type`).
- `String length_field_eval`: INFO field in the evaluation VCF giving each variant's allele length (default `allele_length`).
- `String skip_vep_categories`: Comma-separated VEP consequence categories to skip when comparing annotations (default empty).
- `Int? records_per_shard`: Number of matched variants to keep within a single shard during benchmarking.
- `String? args_string_vcf`: `bcftools view` arguments used to pre-subset the evaluation VCF.
- `String? args_string_vcf_truth`: `bcftools view` arguments used to pre-subset the SNV & indel truth VCF.
- `String? args_string_vcf_sv_truth`: `bcftools view` arguments used to pre-subset the SV truth VCF.
- `String? rename_id_string_vcf`: Expression used to rename variant IDs in the evaluation VCF prior to matching.
- `String? rename_id_string_vcf_truth`: Expression used to rename variant IDs in the SNV & indel truth VCF prior to matching.
- `String? rename_id_string_vcf_sv_truth`: Expression used to rename variant IDs in the SV truth VCF prior to matching.
- `Boolean? rename_id_strip_chr_vcf`: Whether to strip the `chr` prefix when renaming evaluation variant IDs.
- `Boolean? rename_id_strip_chr_vcf_truth`: Whether to strip the `chr` prefix when renaming SNV & indel truth variant IDs.
- `Boolean? rename_id_strip_chr_vcf_sv_truth`: Whether to strip the `chr` prefix when renaming SV truth variant IDs.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `annotations_tsv_benchmark`: TSV mapping evaluation variants to their matched truth variants and match type.
- `benchmark_annotations_summary_tsv`: Per-variant comparison of AF & VEP annotations across matched pairs (only when `compare_annotations` is enabled).
- `benchmark_annotations_stats_tsv`: Summary statistics of the annotation comparison (only when `compare_annotations` is enabled).
- `benchmark_annotations_plots_tarball`: Tarball of annotation-comparison plots (only when `compare_annotations` is enabled).


### [AnnotateCallsetOverlap_AF](wdl/annotation/AnnotateCallsetOverlap_AF.wdl)
TODO


### [AnnotateDbSNP](wdl/annotation/AnnotateDbSNP.wdl)
This workflow annotates each variant in the input VCF with its dbSNP reference SNP identifier (rsID). It matches variants against a per-contig dbSNP VCF on CHROM, POS, REF and ALT, emitting a TSV mapping each matched variant to its `dbSNP_ID`.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `File dbsnp_vcf`: From [references](#references).
- `File dbsnp_vcf_idx`: From [references](#references).

Outputs:
- `annotations_tsv_dbsnp`: TSV mapping variants to their dbSNP identifiers.


### [AnnotateDbVaR](wdl/annotation/AnnotateDbVaR.wdl)
This workflow annotates structural variants in the input VCF with matching records from dbVar. It restricts to variants at or above a minimum length, converts them to a symbolic representation, and matches deletions, duplications and insertions separately against a per-contig dbVar VCF using type-specific size-similarity, reciprocal-overlap and breakpoint-window thresholds. It emits a TSV linking matched variants to their dbVar records.

> **Note:** This workflow considers the DUP region for INS annotated as DUP.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int min_length`: Minimum variant length to consider for matching.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for matching deletions (default `0.7`).
- `Int dup_breakpoint_window`: Breakpoint window, in bp, for matching duplications (default `500`).
- `Float dup_reciprocal_overlap`: Minimum reciprocal overlap for matching duplications (default `0.7`).
- `Float dup_size_similarity`: Minimum size similarity for matching duplications (default `0.7`).
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions (default `0.0`).
- `Float ins_size_similarity`: Minimum size similarity for matching insertions (default `0.5`).
- `File dbvar_vcf`: From [references](#references).
- `File dbvar_vcf_idx`: From [references](#references).

Outputs:
- `annotations_tsv_dbvar`: TSV mapping variants to their matched dbVar records.


### [AnnotateGnomADSTR](wdl/annotation/AnnotateGnomADSTR.wdl)
This workflow annotates tandem-repeat variants in the input VCF with overlapping loci from the gnomAD V4 tandem-repeat catalog. It subsets to tandem-repeat calls and matches each against the catalog using a minimum reciprocal-overlap threshold, emitting a TSV linking calls to their gnomAD TR locus.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Float trv_reciprocal_overlap`: Minimum reciprocal overlap between a call and a catalog locus to be matched (default `0.7`).
- `File gnomad_tr_json`: From [references](#references).

Outputs:
- `annotations_tsv_gnomad_str`: TSV mapping tandem-repeat calls to their gnomAD TR loci.


### [AnnotateGQMetrics](wdl/annotation/AnnotateGQMetrics.wdl)
This workflow computes binned distributions of genotype-quality metrics across the carriers of each variant. For every configured FORMAT field it counts the genotypes whose value falls into each bin — optionally restricted to a variant filter and respecting whether larger or smaller values of that field are better — and can additionally bin allele-balance values. It emits a per-variant TSV of these distribution counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Array[String] gq_fields`: FORMAT fields whose values are binned, one per field to annotate.
- `Array[Array[Int]] gq_bins`: Bin edges for each field in `gq_fields`.
- `Array[String] gq_variant_filters`: Per-field expression restricting which variants the field is binned over.
- `Array[Boolean] gq_larger_field`: Per-field flag indicating whether larger values of the field are better.
- `Boolean ab_annotation`: Whether to additionally compute allele-balance distributions.
- `Array[Float] ab_bins`: Bin edges for allele-balance values.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_gq`: TSV of per-variant genotype-quality (and optional allele-balance) distribution counts.


### [AnnotateIndelTRs](wdl/annotation/AnnotateIndelTRs.wdl)
This workflow flags short insertions and deletions that represent tandem repeats. Using the str-analysis `filter_vcf_to_tandem_repeats` tool, it inspects each indel's sequence and marks it as a tandem repeat when it meets a minimum total repeat length, minimum number of repeats and minimum repeat-unit length, emitting a TSV of the flagged variants.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int min_tandem_repeat_length`: Minimum total tandem-repeat length for an indel to be flagged (default `9`).
- `Int min_repeats`: Minimum number of repeats for an indel to be flagged (default `3`).
- `Int min_repeat_unit_length`: Minimum repeat-unit length for an indel to be flagged (default `1`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `annotations_tsv_trs`: TSV of indels flagged as tandem repeats.


### [AnnotateInSilicoPredictors](wdl/annotation/AnnotateInSilicoPredictors.wdl)
This workflow annotates SNVs and indels with precomputed in-silico predictor scores — CADD, Pangolin, PhyloP, REVEL and SpliceAI — drawn from the gnomAD V4 Hail Tables. It shards the VCF and uses a Hail-based script to look up each variant's scores, emitting a TSV of per-variant predictions.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String annotate_in_silico_predictors_script`: Path to the Hail script that performs the lookups (defaults to the `lr-annotation` repository copy).
- `String genome_build`: Genome build to annotate against (default `GRCh38`).
- `String cadd_ht`: From [references](#references).
- `String pangolin_ht`: From [references](#references).
- `String phylop_ht`: From [references](#references).
- `String revel_ht`: From [references](#references).
- `String spliceai_ht`: From [references](#references).

Outputs:
- `annotations_tsv_insilico`: TSV of per-variant in-silico predictor scores.


### [AnnotateL1MEAID](wdl/annotation/AnnotateL1MEAID.wdl)
This workflow first runs _RepeatMasker_ on the insertions in an input VCF. It then uses its output to run [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) in order to identify, annotate and filter mobile element insertion (MEI) calls. It restricts to insertions at or above a minimum length and emits a TSV of the resulting MEI annotations.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int min_length`: Minimum insertion length to consider for MEI annotation.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_l1meaid`: TSV of L1ME-AID and INTACT_MEI MEI annotations.


### [AnnotateMEDs](wdl/annotation/AnnotateMEDs.wdl)
This workflow annotates mobile element deletions (MEDs) by intersecting the deletions in the input VCF against a catalog of known mobile-element loci. Deletions are extracted to BED form and matched to the catalog using size-similarity, reciprocal-overlap, breakpoint-window and sequence-similarity thresholds, producing a TSV of the deletions identified as MEDs.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for a deletion to match a catalog locus (default `0.7`).
- `Float del_sequence_similarity`: Minimum sequence similarity for a deletion to match a catalog locus (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for a deletion to match a catalog locus (default `0.0`).
- `File mei_catalog`: From [references](#references).

Outputs:
- `annotations_tsv_meds`: TSV of deletions identified as mobile element deletions.


### [AnnotateMEIs](wdl/annotation/AnnotateMEIs.wdl)
This workflow consolidates the mobile element insertion calls produced by the [AnnotateL1MEAID](#annotatel1meaid), [AnnotatePALMER](#annotatepalmer) and [AnnotateSVAN](#annotatesvan) workflows into a single harmonized set. It reconciles the three per-tool annotation TSVs, using the SVAN annotation header for typing, to produce a final TSV of MEI calls.

Inputs:
- `File annotations_tsv_l1meaid`: MEI annotation TSV output by `AnnotateL1MEAID`.
- `File annotations_tsv_palmer`: MEI annotation TSV output by `AnnotatePALMER`.
- `File annotations_tsv_svan`: MEI annotation TSV output by `AnnotateSVAN`.
- `File annotations_header_svan`: Annotation header output by `AnnotateSVAN`, used to type the consolidated fields.
- `Array[String] contigs`: Contigs to annotate.

Outputs:
- `annotations_tsv_meis`: Consolidated TSV of mobile element insertion calls.


### [AnnotatePALMER](wdl/annotation/AnnotatePALMER.wdl)
This workflow leverages [PALMER](https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls. Matching is performed per MEI type using type-specific reciprocal-overlap, size-similarity, sequence-similarity, breakpoint-window and minimum-shared-sample thresholds.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File PALMER_vcf`: VCF of PALMER MEI calls to match against.
- `File PALMER_vcf_idx`: Index for `PALMER_vcf`.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Array[String] mei_types`: MEI types to run on - must be a subset of [`ALU`, `SVA`, `LINE` or `HERVK`].
- `Int min_length`: Minimum insertion length to consider for annotation.
- `File rm_out`: _RepeatMasker_ output for the input VCF's insertions.
- `Int rm_buffer`: Padding, in bp, applied around RepeatMasker annotations when matching.
- `Int ins_breakpoint_window_{ALU,SVA,LINE,HERVK}`: Per-type breakpoint window, in bp, for matching (default `500`).
- `Float ins_reciprocal_overlap_{ALU,SVA,LINE,HERVK}`: Per-type minimum reciprocal overlap for matching (default `0.9`).
- `Float ins_sequence_similarity_{ALU,SVA,LINE,HERVK}`: Per-type minimum sequence similarity for matching (default `0.9`).
- `Float ins_size_similarity_{ALU,SVA,LINE,HERVK}`: Per-type minimum size similarity for matching (default `0.9`).
- `Int ins_min_shared_samples_{ALU,SVA,LINE,HERVK}`: Per-type minimum number of shared samples for matching (default `0`).
- `File ref_fai`: From [references](#references).

Outputs:
- `annotations_tsv_palmer`: TSV of insertions annotated with their PALMER `ME_TYPE`.


### [AnnotateRegion](wdl/annotation/AnnotateRegion.wdl)
This workflow annotates each variant with the genomic region class it falls within — simple repeat (`SR`), segmental duplication (`SD`), RepeatMasker region (`RM`) or unique sequence (`US`) — by intersecting it against the corresponding BED panels. It emits a TSV of per-variant `REGION` assignments.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File simple_repeats_bed`: From [references](#references).
- `File seg_dup_bed`: From [references](#references).
- `File repeat_masked_bed`: From [references](#references).

Outputs:
- `annotations_tsv_region`: TSV of per-variant genomic-region assignments.


### [AnnotateSQMetrics](wdl/annotation/AnnotateSQMetrics.wdl)
This workflow recomputes site-level quality metrics for each variant directly from its genotype-level data. It clears any stale allele-specific INFO fields and recalculates Hardy-Weinberg equilibrium, the inbreeding coefficient, the maximum p(allele balance), and the allele-specific quality approximation, quality-by-depth and variant depth from the per-sample DP, PL and AD fields, emitting a per-variant TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_sq`: TSV of recomputed site-level quality metrics.


### [AnnotateSVAN](wdl/annotation/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It processes insertions and deletions separately, first running Tandem Repeat Finder (TRF) on the inserted or deleted sequence of each SV in the input VCF and then running SVAN over the result, before extracting and aligning the annotations into a single TSV. Before extraction, the `DUP_COORD` field produced by SVAN is reformatted: any `flank_`-prefixed relative coordinates are resolved to absolute genomic positions, preserving the original order of comma-separated values.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File vntr_bed`: From [references](#references).
- `File exons_bed`: From [references](#references).
- `File repeats_bed`: From [references](#references).
- `File mei_fa`: From [references](#references).
- `Array[File] mei_fa_indices`: From [references](#references).
- `File ref_fa`: From [references](#references).
- `Array[File] ref_fa_indices`: From [references](#references).

Outputs:
- `annotations_tsv_svan`: TSV of SVAN annotations.
- `annotations_header_svan`: Header lines describing the SVAN annotation fields.


### [AnnotateSVAnnotate](wdl/annotation/AnnotateSVAnnotate.wdl)
This workflow leverages [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SVs through this workflow, ignoring all SNVs and indels, converting each SV to a symbolic representation before annotating it against coding and noncoding panels and extracting the resulting `PREDICTED_` annotations into a TSV.

> **Note:** This workflow considers the DUP region for INS annotated as DUP.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int min_length`: Minimum length for a variant to be treated as an SV and annotated.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File coding_gtf`: From [references](#references).
- `File noncoding_bed`: From [references](#references).

Outputs:
- `annotations_tsv_svannotate`: TSV of SVAnnotate functional-effect annotations.
- `annotations_header_svannotate`: Header lines describing the SVAnnotate annotation fields.


### [AnnotateVEPHail](wdl/annotation/AnnotateVEPHail.wdl)
This workflow leverages [the Ensembl Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It strips genotypes, scatters the VCF into shards, optionally normalizes and splits multiallelics around the VEP call, and uses Hail in order to run this annotation process in a more efficient and scalable manner before concatenating the per-shard annotations into a single TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCF before annotation.
- `String split_vcf_hail_script`: Path to the Hail script used to scatter the VCF (defaults to the `lr-annotation` repository copy).
- `String vep_annotate_hail_python_script`: Path to the Hail script used to run VEP (defaults to the `lr-annotation` repository copy).
- `String genome_build`: Genome build to annotate against (default `GRCh38`).
- `String vep_json_schema`: Hail type schema describing the structure of VEP's JSON output.
- `String normalize_check_ref`: `bcftools norm` `--check-ref` mode used when normalizing (default `w`).
- `Boolean normalize_vcf`: Whether to normalize and split multiallelics around the VEP call (default `false`).
- `Boolean localize_vcf`: Whether to localize the VCF before annotation (default `true`).
- `Boolean has_index`: Whether the input VCF is indexed (default `true`).
- `Boolean get_chromosome_sizes`: Whether to compute chromosome sizes (default `false`).
- `Boolean split_by_chromosome`: Whether to scatter the VCF by chromosome (default `false`).
- `Boolean split_into_shards`: Whether to scatter the VCF into fixed-size shards (default `false`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).
- `File ref_fa_gz`: bgzipped `ref_fa`, from [references](#references).
- `File ref_fai_gz`: Index for `ref_fa_gz`, from [references](#references).
- `File ref_vep_cache`: From [references](#references).

Outputs:
- `annotations_tsv_vep`: TSV of VEP functional-effect annotations.


### [AnnotateVRS](wdl/annotation/AnnotateVRS.wdl)
This workflow annotates each variant with its GA4GH Variant Representation Specification (VRS) attributes using a seqrepo sequence repository. It runs `vrs-annotate` per contig to add the VRS INFO fields, then extracts them into an annotation TSV of five locating columns (CHROM, POS, REF, ALT, ID) followed by a column for each VRS field.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File seqrepo_tar`: From [references](#references).

Outputs:
- `annotations_tsv_vrs`: TSV of per-variant VRS attributes (`VRS_Allele_IDs`, `VRS_Error`, `VRS_Starts`, `VRS_Ends`, `VRS_States`, `VRS_Lengths`, `VRS_RepeatSubunitLengths`).



## Annotation Utilities
### [AddEndTRs](wdl/annotation_utils/AddEndTRs.wdl)
This utility adds an `END` INFO tag to the tandem-repeat records of a VCF, computed per contig, so that downstream tools correctly interpret the span of each TR call. It outputs the updated VCF.

Inputs:
- `File vcf`: VCF to update.
- `File vcf_idx`: Index for VCF to update.
- `Array[String] contigs`: Contigs to process within the input VCF.

Outputs:
- `vcf_with_end`: VCF with `END` tags added to tandem-repeat records.
- `vcf_with_end_idx`: Index for the updated VCF.


### [AnnotateAlleleType](wdl/annotation_utils/AnnotateAlleleType.wdl)
This utility sets the `allele_type` INFO field on variants in a VCF using three annotation TSVs — one for mobile element deletions, one for mobile element insertions and one for duplications — applying each in turn. Each annotation source can have its values transformed via an optional prefix, suffix and lowercasing. It outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File med_tsv`: TSV of mobile element deletion allele types.
- `File mei_tsv`: TSV of mobile element insertion allele types.
- `File dup_tsv`: TSV of duplication allele types.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `String? med_prefix`: Prefix prepended to mobile element deletion allele-type values.
- `String? med_suffix`: Suffix appended to mobile element deletion allele-type values.
- `Boolean? med_lowercase`: Whether to lowercase mobile element deletion allele-type values.
- `String? mei_prefix`: Prefix prepended to mobile element insertion allele-type values.
- `String? mei_suffix`: Suffix appended to mobile element insertion allele-type values.
- `Boolean? mei_lowercase`: Whether to lowercase mobile element insertion allele-type values.
- `String? dup_prefix`: Prefix prepended to duplication allele-type values.
- `String? dup_suffix`: Suffix appended to duplication allele-type values.
- `Boolean? dup_lowercase`: Whether to lowercase duplication allele-type values.

Outputs:
- `allele_type_annotated_vcf`: VCF annotated with `allele_type`.
- `allele_type_annotated_vcf_idx`: Index for the annotated VCF.


### [AnnotateVcf](wdl/annotation_utils/AnnotateVcf.wdl)
This utility applies a list of annotation TSVs to a VCF as new INFO fields, adding the specified field names, descriptions, types and numbers for each TSV in turn. Each annotation source can optionally have its TSV sorted, the VCF pre-subset and its TSV rows filtered beforehand. It outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[File] annotations_tsvs`: Annotation TSVs to apply, each as a set of INFO fields.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation. When set, each contig VCF is sharded by record count, annotated in parallel and concatenated.
- `Array[Array[String]] info_names`: INFO field names added by each annotation TSV.
- `Array[Array[String]] info_descriptions`: INFO field header descriptions for each annotation TSV.
- `Array[Array[String]] info_types`: INFO field types for each annotation TSV.
- `Array[Array[String]] info_numbers`: INFO field `Number` values for each annotation TSV.
- `Array[Boolean] sort_tsvs`: Per-TSV flag indicating whether to sort the TSV before annotation (default empty).
- `Array[String] subset_vcf_strings`: Per-TSV `bcftools view` arguments used to pre-subset the VCF (default empty).
- `Array[String] awk_tsv_conditions`: Per-TSV `awk` condition used to filter the TSV rows applied (default empty).

Outputs:
- `annotated_vcf`: Annotated VCF.
- `annotated_vcf_idx`: Index for the annotated VCF.


### [AnnotateVcfCleared](wdl/annotation_utils/AnnotateVcfCleared.wdl)
This utility is a variant of [AnnotateVcf](#annotatevcf) that, before applying the annotation TSVs, clears existing annotations and optionally swaps in records from an untrimmed VCF in order to restore full REF/ALT alleles. It then adds the specified INFO fields and outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[File] annotations_tsvs`: Annotation TSVs to apply, each as a set of INFO fields.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Array[Array[String]] info_names`: INFO field names added by each annotation TSV.
- `Array[Array[String]] info_descriptions`: INFO field header descriptions for each annotation TSV.
- `Array[Array[String]] info_types`: INFO field types for each annotation TSV.
- `Array[Array[String]] info_numbers`: INFO field `Number` values for each annotation TSV.
- `File? subset_untrimmed_vcf`: Untrimmed VCF whose records are swapped in to restore full REF/ALT alleles.
- `File? subset_untrimmed_vcf_idx`: Index for `subset_untrimmed_vcf`.
- `Array[Boolean]? sort_tsvs`: Per-TSV flag indicating whether to sort the TSV before annotation.
- `Array[String]? subset_vcf_strings`: Per-TSV `bcftools view` arguments used to pre-subset the VCF.
- `Array[String]? awk_tsv_conditions`: Per-TSV `awk` condition used to filter the TSV rows applied.

Outputs:
- `annotated_vcf`: Annotated VCF.
- `annotated_vcf_idx`: Index for the annotated VCF.


### [CombineCohortTRs](wdl/annotation_utils/CombineCohortTRs.wdl)
This utility combines tandem-repeat VCFs from multiple callers into a single cohort VCF. After checking sample consistency across the inputs, it sets missing filters to pass, tags each caller's calls, assigns TR identifiers, deduplicates overlapping variants and priority-merges the callers on a per-contig basis. It outputs the merged cohort TR VCF.

Inputs:
- `Array[File] tr_vcfs`: Tandem-repeat VCFs to combine, one per caller.
- `Array[File] tr_vcf_idxs`: Indexes for `tr_vcfs`.
- `Array[String] tr_callers`: Caller name for each VCF in `tr_vcfs`, used for tagging and merge priority.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Sample IDs expected across all input VCFs.

Outputs:
- `trgt_merged_vcf`: Merged cohort tandem-repeat VCF.
- `trgt_merged_vcf_idx`: Index for the merged VCF.


### [CombineTRs](wdl/annotation_utils/CombineTRs.wdl)
This utility combines tandem-repeat VCFs from multiple callers for a single sample into one VCF. It checks sample consistency, sets missing filters to pass, tags each caller's calls, assigns TR identifiers, deduplicates overlapping variants and priority-merges the callers per contig. It outputs the combined TR VCF.

Inputs:
- `Array[File] tr_vcfs`: Tandem-repeat VCFs to combine, one per caller.
- `Array[File] tr_vcf_idxs`: Indexes for `tr_vcfs`.
- `Array[String] tr_callers`: Caller name for each VCF in `tr_vcfs`, used for tagging and merge priority.
- `Array[String] contigs`: Contigs to process.
- `String sample_id`: Sample ID expected across all input VCFs.

Outputs:
- `trgt_combined_vcf`: Combined single-sample tandem-repeat VCF.
- `trgt_combined_vcf_idx`: Index for the combined VCF.


### [CombineVcfsAcrossContigs](wdl/annotation_utils/CombineVcfsAcrossContigs.wdl)
This utility concatenates a set of per-contig VCFs into a single VCF, optionally dropping genotypes in the process. It outputs the combined VCF.

Inputs:
- `Array[File] vcfs`: Per-contig VCFs to concatenate.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs corresponding to `vcfs`.
- `Boolean drop_genotypes`: Whether to strip genotypes from the combined VCF (default `false`).

Outputs:
- `concat_vcf`: Combined VCF.
- `concat_vcf_idx`: Index for the combined VCF.


### [CountAnnotations](wdl/annotation_utils/CountAnnotations.wdl)
This utility tallies annotation values across one or more VCFs to produce summary count tables. It always counts at the site level and can optionally count per sample, per allele, per functional gene consequence and as raw value lists, with optional sharding, region splitting and length/region subsetting. It outputs a site-count TSV plus whichever optional breakdowns were requested.

Inputs:
- `Array[File] vcfs`: VCFs whose annotations are counted.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Int? records_per_shard`: Number of variants to keep within a single record-based shard.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding by genomic region instead of record count.
- `Boolean create_per_sample`: Whether to additionally produce per-sample counts (default `false`).
- `Boolean create_per_allele`: Whether to additionally produce per-allele counts (default `false`).
- `Boolean create_list`: Whether to additionally produce raw value-list tables (default `false`).
- `Boolean create_functional`: Whether to additionally produce per-gene functional counts (default `false`).
- `Boolean create_variant_attributes`: Whether to additionally count variant attribute combinations (default `false`).
- `Boolean use_ssd`: Whether to use SSD-backed local disks (default `false`).
- `Boolean split_by_region`: Whether to split each VCF by genomic region before counting (default `false`).
- `String subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCFs (default empty).
- `Int max_length`: Maximum variant length to count, or `-1` for no maximum (default `-1`).
- `Int min_length`: Minimum variant length to count, or `-1` for no minimum (default `-1`).
- `File? ref_fai`: Reference index used when splitting by region.

Outputs:
- `annotation_counts_sites_tsv`: Site-level annotation counts.
- `annotation_counts_functional_tsv`: Per-gene functional counts (when `create_functional`).
- `annotation_counts_functional_samples_tsv`: Per-gene per-sample functional counts (when `create_functional` and `create_per_sample`).
- `annotation_counts_functional_alleles_tsv`: Per-gene per-allele functional counts (when `create_functional` and `create_per_allele`).
- `annotation_counts_list_tsv`: Raw annotation value lists (when `create_list`).
- `annotation_counts_samples_tsv`: Per-sample counts (when `create_per_sample`).
- `annotation_counts_alleles_tsv`: Per-allele counts (when `create_per_allele`).


### [CreateBiallelicVcf](wdl/annotation_utils/CreateBiallelicVcf.wdl)
This utility normalizes a VCF into a streamlined biallelic callset. It splits multiallelic records and left-aligns variants against the reference, sorts the result, adds the `allele_length` and `allele_type` INFO fields, and rewrites each variant ID to `CHROM-POS-REF-ALT` for SNVs or `CHROM-POS-TYPE-LENGTH` otherwise, suffixing any colliding IDs to keep them unique. It outputs the biallelic VCF.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF to process.
- `File ref_fa`: Reference FASTA used for normalization.
- `File ref_fai`: Index for `ref_fa`.

Outputs:
- `biallelic_vcf`: Normalized, sorted biallelic VCF with streamlined variant IDs and `allele_length`/`allele_type` annotations.
- `biallelic_vcf_idx`: Index for the biallelic VCF.


### [CreateCoverageFile](wdl/annotation_utils/CreateCoverageFile.wdl)
This utility builds a binned coverage matrix across a cohort from per-sample mosdepth BED outputs. It tiles the genome into windows, computes the mean coverage and threshold-crossing counts within each bin for every sample, and concatenates the results into a single coverage TSV.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-sample mosdepth coverage BED files.
- `Array[File] mosdepth_bed_indices`: Indexes for `mosdepth_bed_files`.
- `Array[String] contigs`: Contigs over which to compute coverage.
- `Int window_size`: Size, in bp, of each genomic window.
- `Int bin_size`: Size, in bp, of each coverage bin within a window.
- `Array[Int] thresholds`: Coverage thresholds at which to count bins as covered.
- `File ref_fai`: From [references](#references).

Outputs:
- `binned_coverage_tsv`: Binned coverage matrix across the cohort.


### [CreateMetadataFile](wdl/annotation_utils/CreateMetadataFile.wdl)
This utility builds a cohort metadata file by combining a pedigree file with an ancestry-assignment file. It outputs the merged metadata file.

Inputs:
- `File ped_file`: Cohort pedigree file.
- `File ancestry_file`: Two-column file of sample IDs and ancestry labels.

Outputs:
- `metadata`: Merged cohort metadata file.


### [CreatePEDAncestry](wdl/annotation_utils/CreatePEDAncestry.wdl)
This utility generates a minimal pedigree file and an ancestry-assignment file from a list of sample IDs and their sexes. It outputs both files.

Inputs:
- `Array[String] sample_ids`: Sample IDs to include.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.

Outputs:
- `ped`: Generated pedigree file.
- `ancestry`: Generated ancestry-assignment file.


### [CreateReadCountsFile](wdl/annotation_utils/CreateReadCountsFile.wdl)
This utility produces a binned read-counts file for a single sample from its per-contig mosdepth BED outputs, binning counts at a fixed resolution and merging across contigs. It outputs the binned read-counts file.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-contig mosdepth coverage BED files for the sample.
- `Array[File] mosdepth_bed_indices`: Indexes for `mosdepth_bed_files`.
- `Array[String] contigs`: Contigs over which to bin read counts.
- `Int bin_size`: Size, in bp, of each read-count bin.
- `String sample_id`: ID of the sample being processed.
- `File ref_dict`: From [references](#references).

Outputs:
- `binned_read_counts`: Binned read-counts file for the sample.


### [DropGenotypes](wdl/annotation_utils/DropGenotypes.wdl)
This utility strips all genotype (sample) columns from a VCF, optionally sharding by record count for speed. It outputs the resulting sites-only VCF.

Inputs:
- `File vcf`: VCF whose genotypes are dropped.
- `File vcf_idx`: Index for VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `dropped_vcf`: Sites-only VCF.
- `dropped_vcf_idx`: Index for the sites-only VCF.


### [ExtractSampleVcfs](wdl/annotation_utils/ExtractSampleVcfs.wdl)
This utility extracts per-sample VCFs from a cohort VCF, splitting each sample's variants into a SNV/indel VCF and an SV VCF based on a minimum SV length. It outputs the per-sample SNV/indel and SV VCFs.

Inputs:
- `File cohort_vcf`: Cohort VCF to extract from.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[String] contigs`: Contigs to process within the cohort VCF.
- `Array[String] sample_ids`: Samples to extract.
- `Int min_sv_length`: Minimum length at which a variant is routed to the SV VCF rather than the SNV/indel VCF.

Outputs:
- `snv_indel_vcfs`: Per-sample SNV/indel VCFs.
- `snv_indel_vcf_idxs`: Indexes for the SNV/indel VCFs.
- `sv_vcfs`: Per-sample SV VCFs.
- `sv_vcf_idxs`: Indexes for the SV VCFs.


### [FillFormatFields](wdl/annotation_utils/FillFormatFields.wdl)
This utility fills missing FORMAT fields in one VCF using the values from a second, more complete VCF covering the same sites. It supports selectively copying named format fields plus toggles for filling alternate and reference genotypes, unphasing genotypes and adding PL, with optional region sharding. It outputs the refilled VCF.

Inputs:
- `File unfilled_vcf`: VCF whose FORMAT fields are filled.
- `File unfilled_vcf_idx`: Index for `unfilled_vcf`.
- `File filled_vcf`: VCF providing the FORMAT field values.
- `File filled_vcf_idx`: Index for `filled_vcf`.
- `String contig`: Contig to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `Array[String] format_fields`: FORMAT fields to fill.
- `String? include_field`: FORMAT field whose value gates whether a record is filled.
- `String? include_value`: Value of `include_field` required for a record to be filled.
- `Boolean fill_alt_gts`: Whether to fill alternate-allele genotypes (default `false`).
- `Boolean fill_ref_gts`: Whether to fill reference genotypes (default `false`).
- `Boolean unphase_gts`: Whether to unphase genotypes while filling (default `false`).
- `Boolean add_pl`: Whether to add a PL field (default `false`).

Outputs:
- `refilled_vcf`: VCF with FORMAT fields filled.
- `refilled_vcf_idx`: Index for the refilled VCF.


### [FillPhasedGenotypes](wdl/annotation_utils/FillPhasedGenotypes.wdl)
This utility transfers phasing information from a phased VCF onto the genotypes of an unphased VCF over matching sites, optionally sharding each contig by region. It outputs the phased VCF.

Inputs:
- `File phased_vcf`: VCF providing the phasing information.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `File unphased_vcf`: VCF whose genotypes are phased.
- `File unphased_vcf_idx`: Index for `unphased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding each contig.

Outputs:
- `hiphase_phased_vcf`: Phased VCF.
- `hiphase_phased_vcf_idx`: Index for the phased VCF.


### [FilterTRGTCalls](wdl/annotation_utils/FilterTRGTCalls.wdl)
This utility filters a TRGT tandem-repeat VCF, optionally dropping calls below a minimum repeat-unit length or length difference, or above a maximum catalog length. It outputs the filtered VCF.

Inputs:
- `File trgt_vcf`: TRGT VCF to filter.
- `File trgt_vcf_idx`: Index for the TRGT VCF.
- `Int? min_repeat_unit`: Minimum repeat-unit length to retain a call.
- `Int? min_length_diff`: Minimum length difference from the reference to retain a call.
- `Int? max_catalog_length`: Maximum catalog locus length to retain a call.

Outputs:
- `trgt_filtered_vcf`: Filtered TRGT VCF.
- `trgt_filtered_vcf_idx`: Index for the filtered VCF.


### [FindUntrimmedAlleles](wdl/annotation_utils/FindUntrimmedAlleles.wdl)
This utility identifies variants in a VCF whose REF and ALT alleles retain untrimmed shared bases, producing a subset VCF of those records for use in restoring full allele representations downstream. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to scan.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to scan within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during scanning.

Outputs:
- `subset_untrimmed_vcf`: VCF of records with untrimmed alleles.
- `subset_untrimmed_vcf_idx`: Index for the subset VCF.


### [GenerateTRGTJson](wdl/annotation_utils/GenerateTRGTJson.wdl)
This utility converts a TRGT longest-polymer-sequence (LPS) table into per-locus allele-frequency histograms, using a metadata file for sample context, and concatenates them into a single TSV. It outputs the histogram TSV.

Inputs:
- `File lps_tsv`: TRGT LPS table to convert.
- `File metadata_tsv`: Sample metadata used for context.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `trgt_histograms_tsv`: Per-locus allele-frequency histogram TSV.


### [GQCalculateCounts](wdl/annotation_utils/GQCalculateCounts.wdl)
This utility computes GQ-stratified count tables used to derive GQ filtering cutoffs, from both a trio de novo analysis and a truth-set concordance analysis. Counts are bucketed by variant type, allele-length bin and supporting caller. For structural variants (`abs(allele_length) >= 50`), the `CALLER` column expands each call by its supporting callers using the `EV`/`BEV` FORMAT fields written by [SVAddRawCallers](wdl/annotation_utils/SVAddRawCallers.wdl): `kanpig`-backed calls are recorded under `CALLER=kanpig` with their own GQ, calls backed by other callers are split into one row per caller carrying an allelic depth in `EV` (with a per-caller GQ recomputed from that depth), and calls with no `BEV` are recorded with a blank `CALLER`. It outputs one TSV per analysis.

Inputs:
- `Array[File] vcfs`: Cohort VCFs to analyze.
- `Array[File] vcf_idxs`: Indexes for the cohort VCFs.
- `Array[File]? truth_vcfs`: Truth-set VCFs, one per input VCF, for the concordance analysis.
- `Array[File]? truth_vcf_idxs`: Indexes for the truth-set VCFs.
- `Array[Int] length_bins`: Allele-length bin boundaries defining the size buckets.
- `String? subset_vcf_string`: Optional `bcftools view` argument string to pre-subset each VCF.
- `File? ped`: Pedigree used to identify trios for the de novo analysis.
- `File? swap_samples_truth`: Optional sample-swap list applied to the truth VCFs.
- `Boolean run_trio_qc`: Whether to run the trio de novo analysis.
- `Boolean run_truth_qc`: Whether to run the truth-set concordance analysis.
- `Boolean skip_trv`: Whether to skip tandem-repeat variants.
- `Int min_fuzzy_match`: Minimum variant length to perform fuzzy matching for truth concordance (default `20`).
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions during truth concordance (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions during truth concordance (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for matching deletions during truth concordance (default `0.7`).
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions during truth concordance (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions during truth concordance (default `0.0`).
- `Float ins_size_similarity`: Minimum size similarity for matching insertions during truth concordance (default `0.5`).

Outputs:
- `trio_denovo_tsv`: GQ-stratified trio de novo count table.
- `truth_concordance_tsv`: GQ-stratified truth-set concordance count table.


### [IntegrateTRs](wdl/annotation_utils/IntegrateTRs.wdl)
This utility integrates tandem-repeat calls into a base VCF for a cohort. It aligns samples between the base and TR VCFs, sets missing filters to pass, tags TR records with their source catalog, assigns TR identifiers and annotates the base VCF with the integrated TR calls. It outputs the TR-annotated VCF.

Inputs:
- `File vcf`: Base VCF to integrate TRs into.
- `File vcf_idx`: Index for the base VCF.
- `File tr_vcf`: Tandem-repeat VCF to integrate.
- `File tr_vcf_idx`: Index for the TR VCF.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples shared between the base and TR VCFs.
- `Array[File] tr_catalogs`: Catalogs from which the TR calls were derived.
- `Array[String] tr_catalog_ids`: Identifier for each catalog in `tr_catalogs`.

Outputs:
- `tr_annotated_vcf`: Base VCF annotated with integrated TR calls.
- `tr_annotated_vcf_idx`: Index for the annotated VCF.


### [IntegrateVcfs](wdl/annotation_utils/IntegrateVcfs.wdl)
This utility integrates a SNV/indel VCF and an SV VCF into a single cohort VCF. Each input is normalized, harmonized to a common sample set and tagged with a source label and a size-based flag, after which the two are merged and the combined variants are renamed and filtered — for example to flag large SNVs/indels and small SVs. Sample IDs can optionally be swapped first. It outputs the integrated VCF.

Inputs:
- `File snv_indel_vcf`: SNV/indel VCF to integrate.
- `File snv_indel_vcf_idx`: Index for `snv_indel_vcf`.
- `File sv_vcf`: SV VCF to integrate.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples shared between the two VCFs.
- `String snv_indel_vcf_source_tag`: `SOURCE` value applied to the SNV/indel calls.
- `String snv_indel_vcf_size_flag`: Filter flag applied to out-of-range SNV/indel calls.
- `String snv_indel_vcf_size_flag_description`: Header description for `snv_indel_vcf_size_flag`.
- `String sv_vcf_source_tag`: `SOURCE` value applied to the SV calls.
- `String sv_vcf_size_flag`: Filter flag applied to out-of-range SV calls.
- `String sv_vcf_size_flag_description`: Header description for `sv_vcf_size_flag`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `Int min_sv_length`: Length boundary separating SNVs/indels from SVs (default `50`).
- `File? swap_samples_snv_indel`: Sample-ID swap map applied to the SNV/indel VCF.
- `File? swap_samples_sv`: Sample-ID swap map applied to the SV VCF.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `integrated_vcf`: Integrated cohort VCF.
- `integrated_vcf_idx`: Index for the integrated VCF.


### [MergeBackbonePhased](wdl/annotation_utils/MergeBackbonePhased.wdl)
TODO


### [MergeSites](wdl/annotation_utils/MergeSites.wdl)
This utility merges redundant records at the site level within a VCF by collapsing near-identical deletions and insertions. Deletions are collapsed using size-, reciprocal-overlap, sequence- and sample-similarity thresholds plus a breakpoint distance, insertions using size-, sequence- and sample-similarity plus a breakpoint distance, while all other variants pass through untouched. It outputs the merged VCF.

Inputs:
- `File vcf`: VCF to merge.
- `File vcf_idx`: Index for VCF.
- `Int del_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing deletions (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for collapsing deletions (default `0.0`).
- `Float del_sample_similarity`: Minimum sample similarity for collapsing deletions (default `0.5`).
- `Float del_sequence_similarity`: Minimum sequence similarity for collapsing deletions (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for collapsing deletions (default `0.7`).
- `Int del_size_max`: Maximum deletion size to collapse, or `-1` for no maximum (default `-1`).
- `Int del_size_min`: Minimum deletion size to collapse (default `0`).
- `Int ins_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing insertions (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for collapsing insertions (default `0.0`).
- `Float ins_sample_similarity`: Minimum sample similarity for collapsing insertions (default `0.5`).
- `Float ins_sequence_similarity`: Minimum sequence similarity for collapsing insertions (default `0.7`).
- `Float ins_size_similarity`: Minimum size similarity for collapsing insertions (default `0.7`).
- `Int ins_size_max`: Maximum insertion size to collapse, or `-1` for no maximum (default `-1`).
- `Int ins_size_min`: Minimum insertion size to collapse (default `0`).

Outputs:
- `merged_vcf`: Site-merged VCF.
- `merged_vcf_idx`: Index for the merged VCF.


### [MergeVcfs](wdl/annotation_utils/MergeVcfs.wdl)
This utility merges multiple per-contig VCFs covering the same contig into one, handling tandem-repeat and non-tandem-repeat variants separately. Non-TR variants are merged with Truvari using reciprocal-overlap, sequence-, size- and sample-similarity, a breakpoint distance and size bounds, while TR variants are merged on their identifiers, with optional region sharding. It outputs the merged VCF and a merge-summary TSV.

Inputs:
- `Array[File] contig_vcfs`: Per-callset VCFs for the contig being merged.
- `Array[File] contig_vcf_idxs`: Indexes for `contig_vcfs`.
- `String contig`: Contig being merged.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `Int min_truvari_match`: Minimum variant length for Truvari matching (default `20`).
- `Int truvari_breakpoint_window`: Maximum breakpoint distance, in bp, for merging non-TR variants (default `500`).
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for merging non-TR variants (default `0.0`).
- `Float truvari_sample_similarity`: Minimum sample similarity for merging non-TR variants (default `0.0`).
- `Float truvari_sequence_similarity`: Minimum sequence similarity for merging non-TR variants (default `0.7`).
- `Float truvari_size_similarity`: Minimum size similarity for merging non-TR variants (default `0.7`).
- `Int size_min`: Minimum variant size to merge (default `20`).
- `Int size_max`: Maximum variant size to merge (default `50000`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `merged_vcf`: Merged VCF.
- `merged_vcf_idx`: Index for the merged VCF.
- `merge_summary_tsv`: TSV summarizing the merge.


### [ParseAbsoluteOrigin](wdl/annotation_utils/ParseAbsoluteOrigin.wdl)
This utility resolves the relative `ORIGIN` coordinates of duplications and NUMTs into absolute genomic coordinates and annotates them back onto the VCF. `ORIGIN` values prefixed with `flank_` encode coordinates relative to a flanking window and are converted to genome-absolute positions; values already in absolute form are kept as-is. When multiple comma-separated `ORIGIN` values are present — whether flank-relative, absolute, or mixed — each is processed individually and the resulting absolute values are written back in their original order. It outputs the VCF with absolute-origin annotations.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `absolute_origin_vcf`: VCF with absolute `ORIGIN` coordinates.
- `absolute_origin_vcf_idx`: Index for the annotated VCF.


### [PlotPhasingResults](wdl/annotation_utils/PlotPhasingResults.wdl)
This utility evaluates backbone-phasing accuracy by comparing backbone-phased VCFs against base (truth) VCFs. It assigns samples to their base VCFs, compares phased genotypes per contig and aggregates the results into tables broken down by variants outside tandem repeats, TR-enveloped variants and TR variants. It outputs these summary tables plus per-VCF status tables.

Inputs:
- `Array[File] backbone_phased_vcfs`: Backbone-phased VCFs to evaluate.
- `Array[File] backbone_phased_vcf_idxs`: Indexes for `backbone_phased_vcfs`.
- `Array[File] base_vcfs`: Base (truth) VCFs to compare against.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Int max_variants`: Maximum number of variants to evaluate, or `-1` for no limit (default `-1`).
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.

Outputs:
- `outside_tr_table`: Phasing-accuracy table for variants outside tandem repeats.
- `tr_enveloped_table`: Phasing-accuracy table for TR-enveloped variants.
- `trv_table`: Phasing-accuracy table for tandem-repeat variants.
- `missing_samples`: Samples with no matching base VCF.
- `vcf_tables`: Per-VCF variant-status tables.


### [PostProcess](wdl/annotation_utils/PostProcess.wdl)
This utility bundles every genotype-update and post-processing step applied to a near-final callset into one workflow, with a required `run_` Boolean guarding each step so that the input VCF is left untouched when all are set to `false`. The per-record steps are applied in a single pass over the VCF: each variant is first matched against `transfer_vcf` and has its genotypes transferred (when `run_transfer_genotypes` is set) using its unmodified properties, after which the remaining steps — unphasing, ploidy normalization, TR-ID decrementing, MEI pruning, homopolymer flagging, singleton filtering and same-coordinate sorting — run in order. Some steps require an accompanying field — `run_transfer_genotypes` needs `transfer_vcf`, `run_unphase_samples` needs `unphase_samples`, and `run_normalize_ploidy` needs `ped`. The per-record pass can optionally be region-sharded via `shard_bin_size`.

Inputs:
- `File vcf`: VCF to post-process.
- `File vcf_idx`: Index for VCF to post-process.
- `Array[String] contigs`: Contigs to process within the input VCF.
- `Boolean run_transfer_genotypes`: Whether to transfer genotypes from `transfer_vcf` onto heterozygous calls (run first; requires `transfer_vcf`).
- `Boolean run_unphase_samples`: Whether to unphase the samples in `unphase_samples` (requires `unphase_samples`).
- `Boolean run_normalize_ploidy`: Whether to normalize ploidy by sex — clearing chrY female calls, making chrX/chrY male calls hemizygous, enforcing diploidy and right-aligning unphased calls (requires `ped`).
- `Boolean run_decrement_trv_ids`: Whether to decrement tandem-repeat variant IDs.
- `Boolean run_prune_meis`: Whether to reclassify mobile elements whose length falls outside the expected bounds back to plain insertions/deletions.
- `Boolean run_flag_homopolymer_trvs`: Whether to flag tandem repeats with a length-1 shortest motif as `HOMOPOLYMER_TRV`.
- `Boolean run_sorting`: Whether to sort records sharing a coordinate by absolute allele length and variant ID.
- `Boolean run_filter_singletons`: Whether to apply the `SINGLE_READ_SUPPORT` filter to singleton calls.
- `File? transfer_vcf`: VCF whose genotypes are transferred when `run_transfer_genotypes` is set.
- `File? transfer_vcf_idx`: Index for `transfer_vcf`.
- `Array[String] unphase_samples`: Samples to unphase when `run_unphase_samples` is set (defaults to empty).
- `File? ped`: Cohort pedigree file, used for ploidy normalization when `run_normalize_ploidy` is set.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the per-record pass.

Outputs:
- `post_processed_vcf`: Post-processed VCF.
- `post_processed_vcf_idx`: Index for the post-processed VCF.


### [QcAnnotations](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/QcAnnotations.wdl)
This workflow adapts the GATK-SV annotation QC pipeline in order to produce a quality-control report for an annotated callset. It collects VCF-wide site-level statistics per contig, converts the VCF to BED, plots the aggregated site metrics and — when comparison datasets are supplied — benchmarks the callset against them at the site level. It can additionally run a per-sample pass that collects per-sample variant lists, plots per-sample and per-family QC, and benchmarks samples against sample-level comparison datasets, before sanitizing all outputs into a single QC tarball.

Inputs:
- `Array[File] vcfs`: Annotated VCFs to QC.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to QC.
- `File ped_file`: Cohort pedigree file, used for per-family QC.
- `Int variants_per_shard`: Number of variants per shard during VCF-wide collection.
- `Int samples_per_shard`: Number of samples per shard during per-sample collection.
- `Boolean create_variant_attributes`: Whether to compute variant-attribute breakdowns (default `false`).
- `Boolean run_per_sample`: Whether to run the per-sample QC pass (default `true`).
- `Int? subset_sample_count`: Number of samples to subset to before the per-sample pass.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCFs.
- `Int? random_seed`: Random seed for sample subsetting and downsampling.
- `Int? max_gq`: Maximum genotype-quality value used when binning QC metrics.
- `Int? downsample_qc_per_sample`: Number of variants to downsample to per sample during QC.
- `File? sample_renaming_tsv`: TSV mapping sample IDs to renamed IDs prior to QC.
- `Array[Array[String]]? site_level_comparison_datasets`: Site-level datasets to benchmark the callset against.
- `Array[Array[String]]? sample_level_comparison_datasets`: Sample-level datasets to benchmark the callset against.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `sv_vcf_qc_output`: Tarball of the QC plots and metrics.
- `vcf2bed_output`: Merged BED representation of the QC'd VCFs.


### [SplitVcfPerContig](wdl/annotation_utils/SplitVcfPerContig.wdl)
This utility splits a VCF into per-contig VCFs, optionally also producing genotype-free copies and applying fixups such as adding missing INFO header lines, modifying SNV IDs and renaming dbSNP/dbVar contigs. It outputs the per-contig VCFs and their no-genotype counterparts.

Inputs:
- `File vcf`: VCF to split.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to split the VCF into.
- `Boolean create_no_geno`: Whether to also produce genotype-free copies (default `false`).
- `Boolean modify_snv_ids`: Whether to rewrite SNV variant IDs (default `false`).
- `Boolean rename_dbsnp_contigs`: Whether to rename contigs to dbSNP naming (default `false`).
- `Boolean rename_dbvar_contigs`: Whether to rename contigs to dbVar naming (default `false`).
- `Array[String]? missing_info_header_fields`: INFO header lines to add if missing.

Outputs:
- `contig_vcfs`: Per-contig VCFs.
- `contig_vcf_idxs`: Indexes for the per-contig VCFs.
- `contig_no_geno_vcfs`: Per-contig genotype-free VCFs.
- `contig_no_geno_vcf_idxs`: Indexes for the genotype-free VCFs.


### [SubsetTRGTToCatalog](wdl/annotation_utils/SubsetTRGTToCatalog.wdl)
This utility subsets a merged TRGT VCF down to the loci present in a given TRGT catalog BED, per contig. It outputs the catalog-restricted TRGT VCF.

Inputs:
- `File trgt_full_merged_vcf`: Merged TRGT VCF to subset.
- `File trgt_full_merged_vcf_idx`: Index for the TRGT VCF.
- `File trgt_catalog_bed_gz`: bgzipped TRGT catalog BED of loci to retain.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `trgt_merged_vcf`: Catalog-restricted TRGT VCF.
- `trgt_merged_vcf_idx`: Index for the subset VCF.


### [SubsetTsvToColumns](wdl/annotation_utils/SubsetTsvToColumns.wdl)
This utility subsets an annotation TSV to a chosen set of columns, optionally filtering rows to those whose columns match specified values. It outputs the subset TSV.

Inputs:
- `File annotations_tsv`: Annotation TSV to subset.
- `File annotations_header`: Header describing the TSV columns.
- `Array[String] subset_columns`: Columns to retain.
- `Array[Array[String]]? subset_column_values`: Per-column values to filter rows by.

Outputs:
- `subset_tsv`: Column-subset TSV.


### [SubsetVcfToContigs](wdl/annotation_utils/SubsetVcfToContigs.wdl)
This utility subsets a VCF to a chosen set of contigs and concatenates the result. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to subset.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to retain.

Outputs:
- `subset_contigs_vcf`: Contig-subset VCF.
- `subset_contigs_vcf_idx`: Index for the subset VCF.


### [SubsetVcfToPerSample](wdl/annotation_utils/SubsetVcfToPerSample.wdl)
This utility extracts a separate single-sample VCF for each requested sample from a set of cohort VCFs, optionally dropping specified fields first. It outputs the per-sample VCFs.

Inputs:
- `Array[File] cohort_vcfs`: Cohort VCFs to extract from.
- `Array[File] cohort_vcf_idxs`: Indexes for `cohort_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples to extract.
- `String? drop_fields`: Fields to drop from each VCF before extraction.

Outputs:
- `subset_vcfs`: Per-sample VCFs.
- `subset_vcf_idxs`: Indexes for the per-sample VCFs.


### [SVAddRawCallers](wdl/annotation_utils/SVAddRawCallers.wdl)
This utility annotates each SV in a cohort VCF with the set of raw callers that independently support it. For every sample it matches the cohort calls against that sample's per-caller VCFs (Kanpig, cuteSV, Sniffles, Delly, pbsv, Sawfish, dipcall and hapdiff) using reciprocal-overlap, size- and sequence-similarity and a breakpoint window, then merges the support back into the cohort VCF. It outputs the annotated VCF and a TSV of per-caller match counts.

Inputs:
- `File sv_vcf`: Cohort SV VCF to annotate.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] sample_ids`: Samples to process.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `Array[File] kanpig_vcfs`: Per-sample Kanpig VCFs.
- `Array[File] kanpig_vcf_idxs`: Indexes for `kanpig_vcfs`.
- `Array[File?] cutesv_vcfs`: Per-sample cuteSV VCFs.
- `Array[File?] cutesv_vcf_idxs`: Indexes for `cutesv_vcfs`.
- `Array[File?] sniffles_vcfs`: Per-sample Sniffles VCFs.
- `Array[File?] sniffles_vcf_idxs`: Indexes for `sniffles_vcfs`.
- `Array[File?] delly_vcfs`: Per-sample Delly VCFs.
- `Array[File?] delly_vcf_idxs`: Indexes for `delly_vcfs`.
- `Array[File?] pbsv_vcfs`: Per-sample pbsv VCFs.
- `Array[File?] pbsv_vcf_idxs`: Indexes for `pbsv_vcfs`.
- `Array[File?] sawfish_vcfs`: Per-sample Sawfish VCFs.
- `Array[File?] sawfish_vcf_idxs`: Indexes for `sawfish_vcfs`.
- `Array[File?] dipcall_vcfs`: Per-sample dipcall VCFs.
- `Array[File?] dipcall_vcf_idxs`: Indexes for `dipcall_vcfs`.
- `Array[File?] hapdiff_vcfs`: Per-sample hapdiff VCFs.
- `Array[File?] hapdiff_vcf_idxs`: Indexes for `hapdiff_vcfs`.
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `Int truvari_breakpoint_window`: Breakpoint window, in bp, for matching a raw call (default `500`).
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for matching a raw call (default `0.0`).
- `Float truvari_sequence_similarity`: Minimum sequence similarity for matching a raw call (default `0.7`).
- `Float truvari_size_similarity`: Minimum size similarity for matching a raw call (default `0.7`).
- `Int fuzzy_match_breakpoint_window`: Breakpoint window, in bp, for fuzzy-matching a raw call to per-caller stats (default `500`).

Outputs:
- `sv_added_vcf`: Cohort VCF annotated with raw-caller support.
- `sv_added_vcf_idx`: Index for the annotated VCF.
- `sv_match_counts_tsv`: TSV of per-caller match counts.


### [TransformINSToDUP](wdl/annotation_utils/TransformINSToDUP.wdl)
This utility converts insertion variants with `allele_type=dup` into symbolic DUP records where the insertion's duplication source (from `INFO/ORIGIN`) passes two criteria: size similarity between the insertion length and the ORIGIN region length must meet the `dup_size_similarity` threshold, and the insertion POS must fall within the ORIGIN region or within `dup_breakpoint_window` bases of its breakpoints. Passing records have their POS set to the ORIGIN start, REF set to `N`, ALT set to `<DUP>`, and `allele_length` updated to the ORIGIN region length. Variants with multiple comma-separated `ORIGIN` values are skipped and passed through unchanged, as a single unambiguous source region is required for the transformation. Records not meeting the criteria are also passed through unchanged. Supports optional record-count sharding.

Inputs:
- `File vcf`: VCF to transform.
- `File vcf_idx`: Index for `vcf`.
- `Int? records_per_shard`: Number of records per shard for parallel processing.
- `Int dup_breakpoint_window`: Maximum distance (bp) between insertion POS and ORIGIN breakpoints to pass the breakpoint check (default `10`).
- `Float dup_size_similarity`: Minimum size similarity ratio (relative to the larger of the two lengths) between insertion and ORIGIN lengths (default `0.9`).
- `Int min_dup_size`: Minimum insertion size (bp) to consider for transformation (default `50`).

Outputs:
- `transformed_vcf`: VCF with qualifying DUP insertions converted to symbolic DUP records.
- `transformed_vcf_idx`: Index for `transformed_vcf`.


### [VcfToBed](wdl/annotation_utils/VcfToBed.wdl)
This utility converts one or more VCFs to BED format using svtk, optionally including sample and filter columns and splitting BND and complex records, with optional record-count sharding. It outputs the concatenated BED.

Inputs:
- `Array[File] vcfs`: VCFs to convert.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Int? records_per_shard`: Number of variants to keep within a single shard during conversion.
- `Boolean include_samples`: Whether to include the per-sample columns (default `true`).
- `Boolean include_filters`: Whether to include the filter column (default `true`).
- `Boolean split_bnd`: Whether to split BND records (default `false`).
- `Boolean split_cpx`: Whether to split complex records (default `false`).

Outputs:
- `bed`: Concatenated BED representation of the input VCFs.



## Tools
### [Automop](wdl/tools/Automop.wdl)
This tool runs `mop` (via FISS) to clean up unreferenced intermediate files in a Terra workspace, freeing storage. A dry-run mode reports what would be deleted without removing anything.

Inputs:
- `String workspace_namespace`: Terra workspace namespace to clean.
- `String workspace_name`: Terra workspace name to clean.
- `String user`: User running the cleanup.
- `Boolean dry_run`: Whether to report rather than perform deletions.

Outputs:
- `fissfc_log`: Log of the cleanup run.


### [BackbonePhase](wdl/tools/BackbonePhase.wdl)
This tool transfers ("backbone") phasing from a set of base VCFs onto a target VCF for a single contig. It assigns each sample to its base VCF, computes the phase-flip orientation needed to make the target consistent with the backbone, and applies those flips. It outputs the phase-transferred VCF and a list of samples with no matching base VCF.

Inputs:
- `File vcf`: Target VCF to phase.
- `File vcf_idx`: Index for `vcf`.
- `Array[File] base_vcfs`: Base VCFs providing the backbone phasing.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `String contig`: Contig to process.
- `File? swap_samples_base`: Sample-ID swap map applied to the base VCFs.
- `Boolean allow_unphased_match_phase`: Whether to allow unphased genotypes to set the phase orientation (default `false`).

Outputs:
- `transferred_vcf`: Phase-transferred VCF.
- `transferred_vcf_idx`: Index for the transferred VCF.
- `missing_samples`: Samples with no matching base VCF.


### [HiPhase](wdl/tools/HiPhase.wdl)
This tool runs PacBio [HiPhase](https://github.com/PacificBiosciences/HiPhase) to jointly phase a sample's small-variant, SV and (optionally) TRGT VCFs against its aligned reads. It preprocesses and synchronizes the input VCFs per contig, phases them together and optionally haplotags the BAM. It outputs the phased VCF, per-contig phasing statistics and an optional haplotagged BAM.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `File small_vcf`: Small-variant (SNV/indel) VCF to phase.
- `File small_vcf_idx`: Index for `small_vcf`.
- `File sv_vcf`: SV VCF to phase.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] contigs`: Contigs to phase.
- `File? trgt_vcf`: TRGT tandem-repeat VCF to additionally phase.
- `File? trgt_vcf_idx`: Index for `trgt_vcf`.
- `Int? trgt_min_repeat_unit`: Minimum repeat-unit length retained when filtering the TRGT VCF.
- `Boolean? trgt_normalize`: Whether to normalize the TRGT VCF before phasing.
- `Int? trgt_min_length_diff`: Minimum length difference retained when filtering the TRGT VCF.
- `Int? trgt_max_catalog_length`: Maximum catalog length retained when filtering the TRGT VCF.
- `String? hiphase_extra_args`: Additional arguments passed to HiPhase.
- `Boolean run_haplotagging`: Whether to also haplotag the BAM (default `false`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `hiphase_vcf`: Phased VCF.
- `hiphase_vcf_idx`: Index for the phased VCF.
- `hiphase_haplotag_files`: Per-contig haplotag read assignments.
- `hiphase_stats`: Per-contig phasing statistics.
- `hiphase_blocks`: Per-contig phase blocks.
- `hiphase_summary`: Per-contig phasing summaries.
- `hiphase_haplotagged_bam`: Haplotagged BAM (only when `run_haplotagging`).
- `hiphase_haplotagged_bam_idx`: Index for the haplotagged BAM (only when `run_haplotagging`).


### [HiPhaseMerge](wdl/tools/HiPhaseMerge.wdl)
This tool merges per-sample HiPhase-phased VCFs into a cohort VCF on a per-contig basis, optionally also merging the TRGT tandem-repeat calls separately — fixing TRGT `END`/`AL` headers and propagating phase-set tags. It outputs the merged integrated VCF and an optional merged TRGT VCF.

Inputs:
- `Array[File] phased_vcfs`: Per-sample HiPhase-phased VCFs to merge.
- `Array[File] phased_vcf_idxs`: Indexes for `phased_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Boolean merge_trgt`: Whether to additionally merge the TRGT tandem-repeat calls separately.
- `String merge_args`: Arguments controlling the VCF merge (default `--merge id`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `hiphase_merged_integrated_vcf`: Merged integrated cohort VCF.
- `hiphase_merged_integrated_vcf_idx`: Index for the merged integrated VCF.
- `hiphase_merged_trgt_vcf`: Merged TRGT VCF (only when `merge_trgt`).
- `hiphase_merged_trgt_vcf_idx`: Index for the merged TRGT VCF (only when `merge_trgt`).


### [Kanpig](wdl/tools/Kanpig.wdl)
This tool regenotypes a cohort SV VCF against each sample's aligned reads using [Kanpig](https://github.com/ACEnglish/kanpig). It subsets the cohort to the target samples, runs Kanpig per sample with sex-aware ploidy beds, and merges the per-sample genotypes back into both a raw and a processed cohort VCF. It outputs the regenotyped (processed) and raw Kanpig VCFs.

Inputs:
- `File cohort_vcf`: Cohort SV VCF to regenotype.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[File] bams`: Aligned reads, one per sample.
- `Array[File] bais`: Indexes for `bams`.
- `Array[String] sample_ids`: Samples to regenotype.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `String merge_args`: Arguments controlling the per-sample genotype merge (default `--merge id`).
- `String kanpig_params`: Parameters passed to Kanpig (default `--neighdist 500 --gpenalty 0.04 --hapsim 0.97`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).
- `File ploidy_bed_male`: From [references](#references).
- `File ploidy_bed_female`: From [references](#references).

Outputs:
- `sv_kanpig_vcf`: Regenotyped (processed) cohort VCF.
- `sv_kanpig_vcf_idx`: Index for the processed VCF.
- `sv_kanpig_raw_vcf`: Raw Kanpig cohort VCF.
- `sv_kanpig_raw_vcf_idx`: Index for the raw VCF.


### [LRCNVs](wdl/tools/LRCNVs.wdl)
This tool calls copy-number variants across a cohort using the GATK germline CNV (gCNV) pipeline in cohort mode. From per-sample depth profiles over a shared interval list it annotates and filters intervals, determines contig ploidy, fits the gCNV model across scattered interval shards, post-processes the per-sample calls into genotyped interval and segment VCFs, and collects sample- and model-level QC. It outputs the gCNV model, per-sample CNV VCFs, denoised copy ratios and QC status.

Inputs:
- `File intervals`: Interval list over which CNVs are called.
- `Array[String]+ entity_ids`: Sample IDs in the cohort.
- `Array[String]+ depth_profiles`: Per-sample read-depth profiles, aligned to `entity_ids`.
- `String cohort_entity_id`: Identifier for the cohort.
- `Int num_intervals_per_scatter`: Number of intervals processed per scatter shard (default `10000`).
- `File? gatk4_jar_override`: Override GATK4 jar.
- `File? mappability_track_bed`: Mappability track used to annotate intervals.
- `File? mappability_track_bed_idx`: Index for `mappability_track_bed`.
- `File? segmental_duplication_track_bed`: Segmental-duplication track used to annotate intervals.
- `File? segmental_duplication_track_bed_idx`: Index for `segmental_duplication_track_bed`.
- `File? blacklist_intervals`: Intervals to exclude from calling.
- `Int ref_copy_number_autosomal_contigs`: Reference copy number for autosomes (default `2`).
- `Array[String]? allosomal_contigs`: Contigs treated as allosomal.
- `Int maximum_number_events_per_sample`: Maximum number of events permitted per sample (default `1000`).
- The workflow additionally exposes numerous optional gCNV model and contig-ploidy hyperparameters, prefixed `gcnv_` and `ploidy_`, that tune the underlying GATK tasks.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).
- `File ref_dict`: From [references](#references).

Outputs:
- `annotated_intervals`: Intervals annotated with GC content and tracks.
- `filtered_intervals`: Intervals retained after filtering.
- `contig_ploidy_model_tar`: Fitted contig-ploidy model.
- `contig_ploidy_calls_tar`: Per-sample contig-ploidy calls.
- `gcnv_model_tars`: Fitted gCNV models, one per scatter shard.
- `gcnv_calls_tars`: Per-shard per-sample gCNV calls.
- `gcnv_tracking_tars`: Per-shard model-fitting tracking files.
- `genotyped_intervals_vcfs`: Per-sample genotyped interval VCFs.
- `genotyped_segments_vcfs`: Per-sample genotyped segment VCFs.
- `sample_qc_status_files`: Per-sample QC status files.
- `sample_qc_status_strings`: Per-sample QC status strings.
- `model_qc_status_file`: Model-level QC status file.
- `model_qc_string`: Model-level QC status string.
- `denoised_copy_ratios`: Per-sample denoised copy ratios.


### [MethylationProfiling](wdl/tools/MethylationProfiling.wdl)
This tool generates CpG methylation pileups from a haplotagged BAM using [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools), producing combined and per-haplotype methylation BED tracks.

Inputs:
- `File bam`: Haplotagged aligned reads.
- `File bai`: Index for `bam`.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `cpg_combined_bed`: Combined methylation pileup BED.
- `cpg_combined_bed_idx`: Index for the combined BED.
- `cpg_hap1_bed`: Haplotype 1 methylation pileup BED.
- `cpg_hap1_bed_idx`: Index for the haplotype 1 BED.
- `cpg_hap2_bed`: Haplotype 2 methylation pileup BED.
- `cpg_hap2_bed_idx`: Index for the haplotype 2 BED.


### [MinimapAlignment](wdl/tools/MinimapAlignment.wdl)
This workflow leverages [Minimap2](https://github.com/lh3/minimap2) in order to align a sample's maternal and paternal assemblies to a reference, saving the resulting BAMs and PAFs to a specified output directory.

Inputs:
- `File assembly_mat`: Maternal assembly.
- `File assembly_pat`: Paternal assembly.
- `String sample_id`: ID of the sample being aligned.
- `String where_to_save`: Output directory to which the aligned assemblies are saved.
- `String minimap_flags`: Parameters to use when running Minimap2 (default `-a -x asm20 --cs --eqx`).
- `Int minimap_threads`: Number of alignment threads (default `32`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `minimap_assembled_bam_mat`: Aligned maternal-assembly BAM.
- `minimap_assembled_bai_mat`: Index for the maternal BAM.
- `minimap_assembled_paf_mat`: Maternal-assembly PAF alignment.
- `minimap_assembled_bam_pat`: Aligned paternal-assembly BAM.
- `minimap_assembled_bai_pat`: Index for the paternal BAM.
- `minimap_assembled_paf_pat`: Paternal-assembly PAF alignment.


### [MosDepth](wdl/tools/MosDepth.wdl)
This tool runs [mosdepth](https://github.com/brentp/mosdepth) to compute sequencing depth over a sample's BAM per contig, optionally in quantized mode. It outputs the per-base, summary, distribution and (optional) quantized coverage files.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `Array[String] contigs`: Contigs over which to compute depth.
- `Boolean quantize_mode`: Whether to run mosdepth in quantized mode.

Outputs:
- `mosdepth_dist`: Per-contig cumulative coverage distributions.
- `mosdepth_summary`: Per-contig coverage summaries.
- `mosdepth_per_base`: Per-contig per-base coverage.
- `mosdepth_per_base_csi`: Indexes for the per-base coverage.
- `mosdepth_quantized_bed`: Per-contig quantized coverage BEDs (when `quantize_mode`).
- `mosdepth_quantized_bed_csi`: Indexes for the quantized BEDs (when `quantize_mode`).


### [PALMER](wdl/tools/PALMER.wdl)
This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype and then finally integrates these into a final VCF containing multiple MEI types.

Inputs:
- `File bam_pat`: Aligned assembly for paternal haplotype.
- `File bai_pat`: Index for aligned assembly for paternal haplotype.
- `File bam_mat`: Aligned assembly for maternal haplotype.
- `File bai_mat`: Index for aligned assembly for maternal haplotype.
- `Array[String] contigs`: Contigs to run PALMER on.
- `Array[String] mei_types`: Series of MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `String truvari_collapse_params`: Truvari parameters to use when merging across haplotypes.
- `Array[File]? override_palmer_calls_pat`: Optional PALMER calls for paternal haplotype, causing the workflow to bypass its execution.
- `Array[File]? override_palmer_tsd_reads_pat`: Optional PALMER TSD reads for paternal haplotype, causing the workflow to bypass its execution.
- `Array[File]? override_palmer_calls_mat`: Optional PALMER calls for maternal haplotype, causing the workflow to bypass its execution.
- `Array[File]? override_palmer_tsd_reads_mat`: Optional PALMER TSD reads for maternal haplotype, causing the workflow to bypass its execution.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `palmer_pat_calls`: Raw PALMER calls for the paternal haplotype, per MEI type.
- `palmer_pat_tsd_reads`: PALMER TSD reads for the paternal haplotype, per MEI type.
- `palmer_pat_vcfs`: Paternal-haplotype PALMER VCFs, per MEI type.
- `palmer_pat_vcf_idxs`: Indexes for the paternal-haplotype VCFs.
- `palmer_mat_calls`: Raw PALMER calls for the maternal haplotype, per MEI type.
- `palmer_mat_tsd_reads`: PALMER TSD reads for the maternal haplotype, per MEI type.
- `palmer_mat_vcfs`: Maternal-haplotype PALMER VCFs, per MEI type.
- `palmer_mat_vcf_idxs`: Indexes for the maternal-haplotype VCFs.
- `palmer_diploid_vcfs`: Diploid PALMER VCFs merged across haplotypes, per MEI type.
- `palmer_diploid_vcf_idxs`: Indexes for the diploid VCFs.
- `palmer_combined_vcf`: Final VCF combining all MEI types.
- `palmer_combined_vcf_idx`: Index for the combined VCF.


### [PALMERDiploid](wdl/tools/PALMERDiploid.wdl)
This tool runs [PALMER](https://github.com/WeichenZhou/PALMER) on a single sample to generate mobile element insertion calls and convert them to a VCF. It shards the input BAM, runs PALMER per MEI type, merges the shard outputs and converts the raw calls into a per-type VCF, optionally bypassing execution when PALMER calls are supplied directly. It outputs the raw PALMER call and TSD files, per-type VCFs and a combined VCF.

Inputs:
- `File? bam`: Aligned reads to run PALMER on.
- `File? bai`: Index for `bam`.
- `Array[File]? override_palmer_calls`: Optional precomputed PALMER calls, causing the workflow to bypass execution.
- `Array[File]? override_palmer_tsd_files`: Optional precomputed PALMER TSD files, causing the workflow to bypass execution.
- `Array[String] contigs`: Contigs to run PALMER on.
- `String sample`: ID of the sample being processed.
- `String mode`: PALMER run mode.
- `Array[String] mei_types`: MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `palmer_calls`: Raw PALMER calls, per MEI type.
- `palmer_tsd_reads`: PALMER TSD reads, per MEI type.
- `palmer_diploid_vcfs`: Per-MEI-type PALMER VCFs.
- `palmer_diploid_vcf_idxs`: Indexes for the per-type VCFs.
- `palmer_combined_vcf`: Final VCF combining all MEI types.
- `palmer_combined_vcf_idx`: Index for the combined VCF.


### [PALMERMerge](wdl/tools/PALMERMerge.wdl)
This tool merges multiple PALMER MEI VCFs into a single VCF per contig and concatenates the result across contigs. It outputs the merged PALMER VCF.

Inputs:
- `Array[File] vcfs`: PALMER VCFs to merge.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `palmer_merged_vcf`: Merged PALMER VCF.
- `palmer_merged_vcf_idx`: Index for the merged VCF.


### [RepeatMasker](wdl/tools/RepeatMasker.wdl)
This workflow leverages [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated and mobile-element content in the insertions of an input VCF. It extracts each insertion's inserted sequence to a FASTA, optionally restricted to a minimum length, and runs RepeatMasker over it.

Inputs:
- `File vcf`: VCF whose insertions are masked.
- `File vcf_idx`: Index for VCF.
- `Int? min_length`: Minimum insertion length to extract and mask.

Outputs:
- `rm_out`: RepeatMasker output table.
- `rm_fa`: FASTA of the masked insertion sequences.


### [TransferMethylationTags](wdl/tools/TransferMethylationTags.wdl)
This tool transfers methylation base-modification tags (MM/ML) from unaligned BAMs onto an aligned BAM. It extracts the tags per read, then per contig re-attaches them to the aligned reads and sorts, merging the result into a single tagged BAM. It outputs the methylation-tagged BAM and a TSV of the transferred tags.

Inputs:
- `Array[String] unaligned_bam_paths`: Paths to the unaligned BAMs carrying the methylation tags.
- `File aligned_bam`: Aligned BAM to receive the tags.
- `File aligned_bai`: Index for `aligned_bam`.
- `Array[String] contigs`: Contigs to process.
- `Boolean gcs_paths`: Whether `unaligned_bam_paths` are GCS paths (default `false`).
- `String mm_tag`: Base-modification tag name (default `MM`).
- `String ml_tag`: Modification-likelihood tag name (default `ML`).

Outputs:
- `methylation_tagged_bam`: Aligned BAM with methylation tags transferred.
- `methylation_tagged_bai`: Index for the tagged BAM.
- `methylation_tags`: TSV of the transferred methylation tags.


### [TRGT](wdl/tools/TRGT.wdl)
This workflow leverages [TRGT](https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats. 

Inputs:
- `File bam`: Aligned reads.
- `File bai`: Index for aligned reads.
- `String sex`: Sex of sample (one of `M` or `F`).
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).
- `File repeat_catalog_trgt`: From [references](#references).

Outputs:
- `trgt_vcf`: TRGT tandem-repeat genotype VCF.
- `trgt_vcf_idx`: Index for the TRGT VCF.


### [TRGTLPS](wdl/tools/TRGTLPS.wdl)
This tool runs the [trgt-lps](https://github.com/PacificBiosciences/trgt-lps) tool per contig to compute the longest polymer sequence (LPS) within each TRGT-genotyped tandem-repeat locus for every sample, concatenating the results into a single TSV. It outputs the LPS TSV.

Inputs:
- `File vcf`: TRGT VCF to process.
- `File vcf_idx`: Index for the TRGT VCF.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `trgt_lps_tsv`: TSV of per-locus longest polymer sequences.


### [TruvariRemap](wdl/tools/TruvariRemap.wdl)
This tool remaps insertion sequences with minimap2 (via Truvari) in order to flag insertions whose inserted sequence aligns elsewhere in the reference. Each insertion above a minimum length is realigned per contig and assessed against alignment-score and coverage thresholds, emitting a TSV of the remap results.

Inputs:
- `File vcf`: VCF whose insertions are remapped.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to process.
- `Int min_length`: Minimum insertion length to remap.
- `Int max_length`: Maximum insertion length to remap.
- `Int mm2_threshold`: Minimum minimap2 alignment score to flag an insertion.
- `Float cov_threshold`: Minimum alignment coverage to flag an insertion.
- `File ref_fa`: From [references](#references).
- `Array[File] ref_bwa_indices`: BWA indices for `ref_fa`, from [references](#references).

Outputs:
- `annotations_tsv_remap`: TSV of insertion remap results.


### [Vamos](wdl/tools/Vamos.wdl)
This tool runs [Vamos](https://github.com/ChaissonLab/vamos) in order to genotype tandem repeats against a Vamos repeat catalog, in read mode (from an aligned read BAM) and/or assembly mode (from per-haplotype assembly BAMs). It outputs the resulting Vamos VCFs.

Inputs:
- `File? read_bam`: Aligned reads to genotype in read mode.
- `File? read_bai`: Index for `read_bam`.
- `Array[File]? assembly_bams`: Per-haplotype assembly BAMs to genotype in assembly mode.
- `Array[File]? assembly_bais`: Indexes for `assembly_bams`.
- `File repeat_catalog_vamos`: Vamos repeat catalog to genotype against.
- `String sample_id`: ID of the sample being genotyped.

Outputs:
- `vamos_assembly_vcfs`: Per-haplotype assembly-mode Vamos VCFs.
- `vamos_assembly_vcf_idxs`: Indexes for the assembly-mode VCFs.
- `vamos_reads_vcf`: Read-mode Vamos VCF.
- `vamos_reads_vcf_idx`: Index for the read-mode VCF.


### [VcfDist](wdl/tools/VcfDist.wdl)
This tool runs [vcfdist](https://github.com/TimD1/vcfdist) in order to benchmark an evaluation VCF against a truth VCF per contig, computing alignment-based precision/recall and phasing accuracy. It outputs vcfdist's precision-recall, phasing, switch-flip, phase-block and supercluster reports.

Inputs:
- `File vcf_eval`: VCF being evaluated.
- `File vcf_eval_idx`: Index for `vcf_eval`.
- `File vcf_truth`: Truth VCF to evaluate against.
- `File vcf_truth_idx`: Index for `vcf_truth`.
- `Array[String] contigs`: Contigs to evaluate.
- `File? bed_regions`: BED of regions to restrict the evaluation to.
- `String? mode`: vcfdist evaluation mode.
- `Float? threshold`: vcfdist matching threshold.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `File ref_fa`: From [references](#references).

Outputs:
- `vcfdist_phasing_summary_tsv`: Per-contig phasing summaries.
- `vcfdist_switchflips_tsv`: Per-contig switch and flip errors.
- `vcfdist_precision_recall_tsv`: Per-contig precision-recall curves.
- `vcfdist_precision_recall_summary_tsv`: Per-contig precision-recall summaries.
- `vcfdist_phase_blocks_tsv`: Per-contig phase blocks.
- `vcfdist_superclusters_tsv`: Per-contig variant superclusters.
- `vcfdist_query_tsv`: Per-contig query-variant results.
- `vcfdist_truth_tsv`: Per-contig truth-variant results.
- `vcfdist_summary_vcf`: Per-contig annotated summary VCFs.


### [VcfDistCohort](wdl/tools/VcfDistCohort.wdl)
This tool runs [vcfdist](https://github.com/TimD1/vcfdist) across a cohort by pairing each evaluation VCF with its corresponding truth VCF and benchmarking every assigned sample, then aggregating the per-sample results. It outputs cohort-level precision/recall and phasing summaries.

Inputs:
- `Array[File] eval_vcfs`: Evaluation VCFs, one per group.
- `Array[File] eval_vcf_idxs`: Indexes for `eval_vcfs`.
- `Array[File] truth_vcfs`: Truth VCFs, aligned to `eval_vcfs`.
- `Array[File] truth_vcf_idxs`: Indexes for `truth_vcfs`.
- `Array[String] contigs`: Contigs to evaluate.
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `File ref_fa`: From [references](#references).

Outputs:
- `vcfdist_phasing_summary_tsv`: Cohort phasing summary.
- `vcfdist_precision_recall_summary_tsv`: Cohort precision-recall summary.
- `vcfdist_precision_recall_tsv`: Cohort precision-recall curves.
- `vcfdist_switchflips_tsv`: Cohort switch and flip errors.
- `vcfdist_phase_blocks_tsv`: Cohort phase blocks.
- `vcfdist_missing_samples`: Samples with no matching truth VCF.


### [Whatshap](wdl/tools/Whatshap.wdl)
This tool haplotags a sample's BAM against a phased VCF using [WhatsHap](https://github.com/whatshap/whatshap), per contig, then merges the tagged reads into a single BAM. It outputs the haplotagged BAM and per-contig haplotag read lists.

Inputs:
- `File bam`: Aligned reads to haplotag.
- `File bai`: Index for `bam`.
- `File phased_vcf`: Phased VCF used to assign haplotypes.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `String? extra_args`: Additional arguments passed to WhatsHap.
- `File ref_fa`: From [references](#references).
- `File ref_fai`: From [references](#references).

Outputs:
- `haplotagged_bam`: Haplotagged BAM.
- `haplotagged_bai`: Index for the haplotagged BAM.
- `haplotag_lists`: Per-contig haplotag read assignments.



## Output Schema
### INFO Fields
- `allele_length`: Allele length - positive for insertions, negative for deletions and 0 for SNVs.
- `allele_type`: Allele type, which is one of the below.
	- `snv`: Single nucleotide variant.
	- `ins`: Insertion.
	- `del`: Deletion.
	- `trv`: Tandem repeat.
	- `dup`: Tandem duplication.
	- `dup_interspersed`: Interspersed duplication.
	- `complex_dup`: Complex duplication.
	- `inv_dup`: Inverted duplication.
	- `numt`: Nuclear-mitochondrial segment.
	- `{ME_TYPE}_ins`: Mobile element insertion, where `{ME_TYPE}` is one of `ALU`, `LINE` or `SVA`.
	- `{ME_TYPE}_del`: Mobile element deletion, where `{ME_TYPE}` is one of `ALU`, `LINE` or `SVA`.
- Allele Frequencies.
	- `AC`: Count of non-reference alleles.
	- `AF`: Proportion of alleles that are non-reference.
	- `AN`: Count of alleles genotyped.
	- `AP_allele`: Allele purity per-allele (multiallelic sites only).
	- `freq_het`: Heterozygous genotype frequency (biallelic sites only).
	- `freq_homalt`: Homozygous alternate genotype frequency (biallelic sites only).
	- `freq_homref`: Homozygous reference genotype frequency (biallelic sites only).
	- `grpmax`: Population label with maximum allele frequency (biallelic sites only).
	- `LPS_allele`: Longest polymer sequence per-allele (multiallelic sites only).
	- `MC_allele`: Motif count per-allele (multiallelic sites only).
	- `n_bi_genos`: Total number of samples with complete genotypes (biallelic sites only).
	- `NCR`: Proportion of alleles that don't have a genotype call.
	- `nhet`: Number of samples with heterozygous genotypes (biallelic sites only).
	- `nhomalt`: Number of samples with homozygous alternate genotypes (biallelic sites only).
	- `nhomref`: Number of samples with homozygous reference genotypes (biallelic sites only).
- `dbSNP_ID`: Variant ID from dbSNP for matched variants.
- `dbVaR_ID`: Variant ID from dbVaR for matched variants.
- Genotype Quality Metrics.
	- `ab_hist_alt_bin_freq`: Histogram for AB in heterozygous individuals calculated on high quality genotypes; bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00.
	- `dp_hist_all_bin_freq`: Histogram for DP calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
	- `dp_hist_all_n_larger`: Count of DP values falling above the highest histogram bin edge, calculated on high quality genotypes.
	- `dp_hist_alt_bin_freq`: Histogram for DP in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
	- `dp_hist_alt_n_larger`: Count of DP values falling above the highest histogram bin edge in heterozygous individuals, calculated on high quality genotypes.
	- `gq_hist_all_bin_freq`: Histogram for GQ calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
	- `gq_hist_alt_bin_freq`: Histogram for GQ in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
	- `sd_hist_all_bin_freq`: Histogram for SD calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
	- `sd_hist_alt_bin_freq`: Histogram for SD in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
- gnomAD Overlap.
	- `gnomAD_STR`: gnomAD STR ID overlapping this tandem repeat variant.
	- `gnomAD_V4_match_filter`: FILTER associated with matched gnomAD V4 variant.
	- `gnomAD_V4_match_ID`: Variant ID of matched gnomAD V4 variant.
	- `gnomAD_V4_match_source`: Source of matched gnomAD V4 variant, which is one of the below.
		- `SNV_indel`: SNV & indel callset.
		- `SV`: SV callset.
	- `gnomAD_V4_match_type`: Method for generating match with gnomAD V4, which is one of the below.
		- `BEDTOOLS_CLOSEST`: Bedtools closest match finetuned for SVs.
		- `EXACT_MATCH`: Exact match across CHROM, POS, REF and ALT.
		- `TRUVARI_{X}`: Truvari match requiring X% sequence similarity.
- `HOMOPOLYMER_TRV`: Flag indicating a tandem repeat call where the shortest motif has length 1.
- Insilico Predictors.
	- `cadd_phred`: CADD Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger values are more deleterious.
	- `cadd_raw_score`: Raw CADD score, interpretable as the extent to which the annotation profile for a given variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' (positive values). Larger values are more deleterious.
	- `pangolin_largest`: Pangolin's largest delta score across 2 splicing consequences, reflecting the probability of the variant being splice-altering.
	- `phylop`: Base-wise conservation score across the 241 placental mammals in the Zoonomia project. Score ranges from -20 to 9.28, where negative values reflect acceleration (faster evolution than expected under neutral drift) and positive values reflect conservation (slower than expected evolution).
	- `revel_max`: Maximum REVEL score at a site's MANE Select or canonical transcript. An ensemble score for predicting the pathogenicity of missense variants based on 13 other variant predictors. Scores range from 0 to 1, with higher scores predicted to be more likely to be deleterious.
	- `spliceai_ds_max`: Illumina's SpliceAI max delta score, interpreted as the probability of the variant being splice-altering.
- `ORIGIN`: Origin of duplicated sequence for duplications and NUMTs.
- `PREDICTED_*`: Annotations from SVAnnotate, which are all prefixed by `PREDICTED_`.
- `REGION`: Genomic region, which is one of `SR` (for simple repeats), `SD` (for segmental duplications), `RM` (for RepeatMasker annotated regions) or `US` (for unique sequences, or more simply, none of the previous regions).
- Site Quality Metrics.
	- `AS_pab_max`: Allele-specific maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5.
	- `AS_QD`: Allele-specific variant call confidence normalized by depth of sample reads supporting a variant.
	- `AS_QUALapprox`: Allele-specific sum of PL[0] values; used to approximate the QUAL score.
	- `AS_VarDP`: Allele-specific depth over variant genotypes (does not include depth of reference samples).
	- `HWE`: Hardy-Weinberg equilibrium p-value.
	- `inbreeding_coeff`: Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium).
- `SOURCE`: Source of call, which is one of the below.
	- `DeepVariant`: SNV or indel call made by the DeepVariant pipeline.
	- `HPRC_SV_Integration`: Structural variant call made by the HPRC SV Integration pipeline.
	- `TRExplorer`: Tandem repeat loci from the TRExplorer v1.0.1 catalog.
	- `Vamos`: Tandem repeat loci from the Vamos v2.1 catalog.
- `SUB_FAMILY`: Sub-family for MEI calls.
- Tandem Repeats.
	- `MOTIFS`: Motifs that the tandem repeat is composed of.
	- `STRUC`: Structure of the tandem repeat region.
	- `TRID`: TR identifier for TR calls; source of enveloping variant with `allele_type = "trv"` for non-TR calls with the `TR_ENVELOPED` flag.
- `TR_ENVELOPED`: Flag indicating a variant with `allele_type != "trv"` is completely enveloped by a variant with `allele_type = "trv"`.
- `TR_PARSED`: Flag indicating a variant with `allele_type != "trv"` is flagged as a tandem repeat by the `AnnotateIndelTRs` workflow.
- `vep`: Annotations from the Variant Effect Predictor (VEP).
- `VRS_Allele_IDs`: Computed identifiers for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles [VRS version=2.0.1;VRS-Python version=2.3.1].
- `VRS_Ends`: Interresidue coordinates used as the location ends for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles.
- `VRS_Error`: Error message if an error occurred computing a VRS Identifier.
- `VRS_Lengths`: Length values from ReferenceLengthExpression states for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles.
- `VRS_RepeatSubunitLengths`: Repeat subunit length values from ReferenceLengthExpression states for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles.
- `VRS_Starts`: Interresidue coordinates used as the location starts for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles.
- `VRS_States`: Literal sequence states used for the GA4GH VRS Alleles corresponding to the GT indexes of the REF and ALT alleles.

### FILTER Values
- `LARGE_SNV_INDEL`: Variant with `SOURCE = "DeepVariant"` that has  `INFO/allele_length ≥ 50`.
- `MONOALLELIC`: Site represents one ALT allele in a region with multiple variants that could not be unified into non-overlapping multi-allelic sites.
- `SINGLE_READ_SUPPORT`: Variant supported by a single read in a single sample.
- `SMALL_SV`: Variant with `SOURCE = "HPRC_SV_Integration"` that has `INFO/allele_length < 50`.



## Code Conventions
### WDL
- Workflows should be structured in the following order, with each of the below separated by a blank line:
	1. Imports.
	2. Inputs.
	3. Definition of variables dynamically generated in the workflow itself.
	4. Calls to tasks.
	5. Outputs.
- Tasks should be structured in the following order, with each of the below separated by a blank line:
	1. Inputs.
	2. Definition of variables dynamically generated in the task itself.
	3. Command.
	4. Outputs.
	5. Runtime settings - default parameters, followed by a select first with the runtime override, then the actual runtime block.
- Inputs should be structured in the following order, with each of the below separated by a blank lines:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable). Also the contigs to be run on as well as the prefix.
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. Reference files - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime-related information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
	5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating outputs.
- Workflow imports should not be renamed using the `as` operator.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain be any consecutive blank lines - i.e. they should have a maximum of one blank line at a time.
- Inputs passed to a task should not have blank lines between inputs.
- The order of inputs passed to a task should reflect their order in the inputs on the workflow level.
- Inputs passed to a task should have a space on either side of the `=` character.
- The inputs section of a task should not have blank lines between inputs.
- Tasks should always have input fields `docker` and `runtime_attr_override` defined, though what is passed to each one of these when calling the task should be explicitly named - e.g. `gatk_docker` and `runtime_attr_override_svannotate` respectively.
- Tasks should also have a prefix input defined, which is passed and set at the workflow level - the outputs from the task should simply use the prefix along with the file type.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_size` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_size`.
- The default `mem_gb`, `boot_disk_gb` and `compute_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 2, besides tasks in `Helpers.wdl` which should be 1.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- The names of inputs, variables and outputs should include a `_` to separate words, and be entirely lowercase unless they refer to a noun that is capitalized (e.g. PALMER or L1MEAID) - i.e. they should always be in snake case.
- There should never be any additional indentation in order to better align parts of the code to the length or horizontal/vertical spacing of other components in its section - indentation should only be applied at the start of a line.
- All mentions of `fasta` should instead use `fa` - e.g. `ref_fa` instead of `ref_fasta`.
- All mentions of `fasta_index`, `fasta_fai` or `fa_fai`  should instead use `fai` - e.g. `ref_fai` instead of `ref_fasta_index`, `ref_fasta_fai` or `ref_fa_fai`.
- All mentions of `vcf_index` or `vcf_tbi` should instead use `vcf_idx`.
- All VCFs should have suffix `_vcf`, and be coupled with a VCF index file that has a suffix `_vcf_idx`.
- Tasks that can be generalized and used across workflows should live in `Helpers.wdl` and be imported by consumer workflows, rather than explicitly defined in a standalone workflow itself.
- Workflow file names must always match the workflow defined within them.
- Annotation workflows should always output a TSV file rather than a VCF, unless its annotations are done for every single variant in the input VCF or if the underlying workflow is designed to annotate variants in a VCF. 


### Python
- All code should be formatted in-line with black's formatting, which can be applied via `black .`.
- All code should be compliant with `flake8`.


### Codebase
- Workflows in `wdl/annotation/` should begin with _Annotate_.
- Workflows directly run in the pipeline should be in one of `wdl/annotation/`, `wdl/annotation_utils/` or `wdl/tools/`, and should have entries in `dockstore.yml` and `README.md`.


### Workspace
- All reference files - i.e. those not specific to an input callset - should be passed in via workspace data.
- All dockers should be passed in via workspace data.



## Processing Notes
Callset Regeneration.
- _FillFormatFields_ on _allele_type_annotated_vcf_.
- _ParseAbsoluteOrigin_.
- Annotation: _AnnotateAgeMetrics_, _AnnotateCallsetOverlap_, _AnnotateDbVaR_, _AnnotateGQMetrics_, _AnnotateSQMetrics_, _AnnotateSVAnnotate_
- _AnnotateVcf_Functional_.
- _AnnotateAF_.
- _AnnotateVcf_Downstream_.
- _PostProcess_.
- Untrim variants: 
	- _FindUntrimmedAlleles_.
	- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbSNP_, _AnnotateDbVaR_, _AnnotateInSilicoPredictors_, _AnnotateVRS_.
	- _AnnotateVcfCleared_.
