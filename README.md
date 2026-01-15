# Long-Read Annotation
This repository serves as a home for all scripts, workflows and processes for annotating long-read callsets.



## Cohort
- HGSVC3.
	- Data:
		- 65 total samples.
		- 32x aligned reads, produced after downsampling by Fabio.
		- High coverage assemblies, derived directly from HGSVC.
		- [VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/).
		- [Samples](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
	- Notes:
		- Metadata file has 67 samples, though it misses NA24385 (HG002 in Terra) from its VCF, renames NA21487 from its VCF to GM21487 and additionally includes GM19320, GM20355 & GM19129.
- HPRC2.
	- Data:
		- 232 total samples.
		- High coverage aligned reads, produced directly by Fabio.
		- High coverage assemblies, derived directly from HPRC.
		- Base VCFs: From Jiadong & HPRC teams.
		- VCF: Produced from base VCFs using _IntegrateVcfs_.
		- [Samples](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
	- Notes: 
		- Metadata file has 234 samples, as it also includes GRC38 and CHM13.
		- HPRC's base VCF has 233 samples, including an additional sample HG03492.
		- Jiadong's base VCF has 1218 samples, including many additional samples.
- Additional Notes:
	- Overlapping Samples: 5 (HG002, HG00733, HG02818, NA19036, NA19240).
	- Mismatched Samples: NA24385 in HGSVC3 is named HG002 in HPRC.
	- The pedigrees in the Terra workspace were made to match the callset VCFs.



## References
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.
- `contigs`: List of chr1 to chr22, plus chrX and chrY.
- `contigs_primary`: List of chr1 to chr22.
- `exons_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `mei_fa`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `mei_fa_amb`: Index for `mei_fa`.
- `mei_fa_ann`: Index for `mei_fa`.
- `mei_fa_bwt`: Index for `mei_fa`.
- `mei_fa_mmi`: Index for `mei_fa`.
- `mei_fa_pac`: Index for `mei_fa`.
- `mei_fa_sa`: Index for `mei_fa`.
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the GATK-SV featured workspace.
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the GATK-SV featured workspace.
- `ref_fa`: Reference sequence, derived from what was used for [PAV](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/).
- `ref_fai`: Index for `ref_fa`.
- `ref_vep_cache`: [v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/). This also contains a series of additional references, which include the MANE protein coding GTF, GENCODE gene list, ClinVar annotation set etc, with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).
- `repeats_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `repeat_catalog_trgt`: [Panel for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [Ben's repository](https://github.com/broadinstitute/tandem-repeat-catalog/releases) as used in All of Us Phase 2.
- `top_level_fa`: [hg38 Release 76](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GRCh38.dna.toplevel.chr.fa.gz) from [Ensemble](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `vntr_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.



## Dockers
- `hail_docker`: `us.gcr.io/talkowski-sv-gnomad/shineren:hail` from Eren.
- `gatk_docker`: `us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT` from the GATK-SV Featured Workspace.
- `minimap_docker`: `eichlerlab/assembly_eval:0.2` from the Eichler lab.
- `minimap_finalize_docker`: `us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3` from the DSP long-reads team.
- `repeatmasker_docker`: `dfam/tetools:1.8` from DFam.
- `sv_base_mini_docker`: `us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52` from the GATK-SV Featured Workspace.
- `trgt_docker`: `us.gcr.io/broad-dsp-lrma/lr-trgt:3.0.0` from the DSP long-reads team.



## Primary Annotation Workflows
### [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow leverages [AnnotateVcf](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling/workflows/broad-firecloud-dsde-methods/20-AnnotateVcf) from the GATK-SV pipeline in order to annotate internal allele frequencies based on sample sexes and ancestries. It runs on all variants in the input VCF, including SVs.

Note that this workflow directly annotates the input VCF rather than outputting a TSV of annotations.

Inputs:
- `sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `ped_file`: Six column file containing the cohort pedigree, with specifications described in [this article](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
- `contigs`.
- `par_bed`.


### [AnnotateL1MEAID](wdl/annotation/AnnotateL1MEAIDFilter.wdl)
This workflow first runs _RepeatMasker_ on an input VCF. It then uses its output to run [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) in order to annotate and filter MEI calls.

TODO


### [AnnotateMEDs](wdl/annotation/AnnotateMEDs.wdl)
TODO


### [AnnotatePALMER](wdl/annotation/AnnotatePALMER.wdl)
This workflow leverages [PALMER](https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls.

Inputs:
- `rm_fa`: Output by _RepeatMasker_.
- `rm_out`: Output by _RepeatMasker_.
- `contigs`.
- `ref_fai`.


### [AnnotateSVAN](wdl/annotation/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

Inputs:
- `contigs`.
- `exons_bed`.
- `mei_fasta`.
- `ref_fa`.
- `repeats_bed`.
- `vntr_bed`.


### [AnnotateSVAnnotate](wdl/annotation/AnnotateSVAnnotate.wdl)
This workflow leverages [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SV through this workflow, ignoring all SNVs and InDels.

Note that this workflow directly annotates the input VCF rather than outputting a TSV of annotations.

Inputs:
- `coding_gtf`.
- `contigs`.
- `noncoding_bed`.


### [AnnotateVEPHail](wdl/annotation/AnnotateVEPHail.wdl)
This workflow leverages [the Ensembl Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It requires numerous references that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

Inputs:
- `ref_fa`.
- `ref_vep_cache`.
- `top_level_fa`.


### [BenchmarkAnnotations](wdl/annotation/BenchmarkAnnotations.wdl)
This workflow ingests two VCFs and finds matching variants across them in order to compare the AF & VEP annotations of these matched pairs. This serves as a degree of benchmarking, as it ensures that annotations applied to a larger cohort (e.g. gnomAD) are in line with those we annotate. It also enables the identification of variants that are outliers relative to exiting cohorts by pulling out those with a large amount of discordance in their annotation across the callsets.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs:
1. Exact match across CHROM, POS, REF and ALT.
2. Truvari match with overlap percentages of 90%, 70% and 50%.
3. Matching based on `bedtools closest`, finetuned for SVs.

Inputs:
- `vcf_truth`: VCF containing SNV & indels to evaluate against.
- `vcf_sv_truth`: VCF containing SVs to evaluate against.
- `ref_fa`.
- `ref_fai`.



## Downstream Annotation Workflows
### [AnnotateVcf](wdl/annotation_utils/AnnotateVcf.wdl)
TODO


### [IntegrateVcfs](wdl/annotation_utils/IntegrateVcfs.wdl)
TODO


### [MergeVEPAF](wdl/annotation_utils/MergeVEPAF.wdl)
TODO



## Additional Tools
### [MinimapAlignment](wdl/tools/MinimapAlignment.wdl)
This workflow leverages [Minimap2](https://github.com/BeckLaboratory/agglovar) in order to align assemblies to a reference.

Inputs:
- `assembly_mat`: Maternal assembly.
- `assembly_pat`: Paternal assembly.
- `minimap_flags`: Parameters to use when running Minimap2.
- `ref_fa`.
- `ref_fai`.

### [PALMER](wdl/tools/PALMER.wdl)
This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype and then finally integrates these into a final VCF containing multiple MEI types.

Inputs:
- `bam_pat`: Aligned assembly for paternal haplotype.
- `bai_pat`: Index for aligned assembly for paternal haplotype.
- `bam_mat`: Aligned assembly for maternal haplotype.
- `bai_mat`: Index for aligned assembly for maternal haplotype.
- `mei_types`: Series of MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `truvari_collapse_params`: Truvari parameters to use when merging across haplotypes.
- `override_palmer_calls_pat`: Optional PALMER calls for paternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_tsd_reads_pat`: Optional PALMER TSD reads for paternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_calls_mat`: Optional PALMER calls for maternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_tsd_reads_mat`: Optional PALMER TSD reads for maternal haplotype, causing the workflow to bypass its execution. 
- `contigs`.
- `ref_fa`.
- `ref_fai`.


### [PALMERDiploid](wdl/tools/PALMERDiploid.wdl)
TODO


### [PALMERMerge](wdl/tools/PALMERMerge.wdl)
TODO


### [RepeatMasker](wdl/tools/RepeatMasker.wdl)
This workflow leverages [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated regions in an input VCF.


### [TRGT](wdl/tools/TRGT.wdl)
This workflow leverages [TRGT](https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats. 

Inputs:
- `bam`: Aligned reads.
- `bai`: Index for aligned reads.
- `sex`: Sex of sample (one of `M` or `F`).
- `ref_fa`.
- `ref_fai`.
- `repeat_catalog_trgt`.


### [TRGTMerge](wdl/tools/TRGTMerge.wdl)
TODO



## Code Conventions
### WDL
- Workflows should be structured in the following order, with each of the below separated by a blank line:
	1. Imports.
	2. Inputs.
	3. Definition of variables dynamically generated in the workflow itself - unless these require outputs from other tasks, in which case they can go after the tasks they are dependent on.
	4. Calls to tasks.
	5. Outputs.
- Tasks should be structured in the following order, with each of the below separated by a blank line:
	1. Inputs.
	2. Definition of variables dynamically generated in the task itself.
	3. Command.
	4. Outputs.
	5. Runtime settings - which should first be its default runtime settings, followed by a select first with the runtime override, then the actual runtime block.
- Workflow inputs should be structured in the following order, with each of the below separated by a blank lines:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable).
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. References - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
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
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_size` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_size`.
- The default `mem_gb`, `boot_disk_gb` and `compute_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 1.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- The names of inputs, variables and outputs should include a `_` to separate words, and be entirely lowercase unless they refer to a noun that is capitalized (e.g. PALMER or L1MEAID) - i.e. they should always be in snake case.
- There should never be any additional indentation in order to better align parts of the code to the length or horizontal/vertical spacing of other components in its section - indentation should only be applied at the start of a line.
- All mentions of `fasta` should instead use `fa` - e.g. `ref_fa` instead of `ref_fasta`.
- All mentions of `fasta_index`, `fasta_fai` or `fa_fai`  should instead use `fai` - e.g. `ref_fai` instead of `ref_fasta_index`, `ref_fasta_fai` or `ref_fa_fai`.
- All mentions of `vcf_index` or `vcf_tbi` should instead use `vcf_idx`.
- All VCFs should have suffix `_vcf`, and be coupled with a VCF index file that has a suffix `_vcf_idx`.
- Tsks that can be generalized and used across workflows should live in `Helpers.wdl` and be imported by consumer workflows, rather than explicitly defined in a standalone workflow itself.
- Workflow file names must always match the workflow defined within them.
- All code written should use 4-space tabs for indentation.


### Python
- All code should be formatted in-line with black's formatting, which can be applied via `black .`.
- Running `flake8` should yield no errors.
- All code written should use 2-space tabs for indentation.


### Workspace Setup
- All references should be passed in via workspace data.
- All dockers should be passed in via workspace data.
