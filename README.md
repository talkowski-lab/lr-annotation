# Long-Read Annotation
Tools for callset annotation in long-read cohorts.


##  Short Tandem Repeats (STRs)
The [AnnotateSTRs.wdl](wdl/AnnotateSTRs.wdl) workflow is based on a [script](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/filter_vcf_to_STR_variants.py) developed by Ben Weisburd to annotate STRs based on their sequence context. It involves reading the nucleotide content of each variant, and comparing this to its surrounding reference genome context.

References:
- `STR Analysis Package`: [Current version](https://github.com/broadinstitute/str-analysis/tree/main) from the Github repository, which is built directly in the [docker](dockerfiles/Dockerfile.AnnotateSTRs).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.


## Insertions
The [AnnotateSVAN.wdl](wdl/AnnotateSVAN.wdl) workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

References:
- `SVAN Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly in the [docker](dockerfiles/Dockerfile.AnnotateSVAN).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `vntr_bed`: From the recommended [hg38 inputs](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the repository.
- `exons_bed`: From the recommended [hg38 inputs](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the repository.
- `repeats_bed`: From the recommended [hg38 inputs](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the repository.
- `mei_fasta`: From the recommended [hg38 inputs](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the repository.


## Functional Consequences
### Small Variants (< 50 bp)
The [AnnotateVEPHail.wdl](wdl/AnnotateVEPHail.wdl) leverages the [Ensembl Variant Effect Predictor (Ensembl VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) to annotate predicted functional effects based on site-level information. It requires numerous reference files that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

References:
- `VEP Package`: [v105](https://github.com/REPBIO-LAB/SVAN) from the VEP Dockerhub repository, which is used as the base image in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `Hail Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly using `pip` in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `top_level_fa`: [Ensemble hg38 - Release 76](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `alpha_missense_file`: [Alpha Missense hg38](https://zenodo.org/records/8208688).
- `eve_data`: [EVE hg38](https://evemodel.org/api/proteins/bulk/download/).
- `ref_vep_cache`: [v105](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).
- Additional References: Based on the pre-defined references used in the installed version of VEP, with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). These additional references include the MANE protein coding GTF, GENCODE gene list, ClinVar annotation set etc.


### Structural Variants (≥ 50 bp)
The [AnnotateSVAnnotate.wdl](wdl/AnnotateSVAnnotate.wdl) workflow employs the GATK tool [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) to annotate predicted functional effects for structural variants. It conditionally only runs structural variants through this workflow, ignoring all SNVs and InDels.

References:
- `GATK Package`: [v4.5.0.0](https://github.com/broadinstitute/gatk), which is built in the base image `quay.io/ymostovoy/lr-utils-basic:2.0` that is used in the [docker](dockerfiles/Dockerfile.AnnotateSVAnnotate).
- `coding_gtf`: [Gencode v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from gnomAD.
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the GATK-SV featured workspace.


## Internal Allele Frequency (AF)
The [AnnotateVcf.wdl](https://github.com/broadinstitute/gatk-sv/blob/kj_gnomad_lr/wdl/AnnotateVcf.wdl) workflow, which is located in the `kj_gnomad_lr` branch of the GATK-SV repository, annotates the internal allele frequencies based on a passed list of sample annotations. It runs all variants in the input VCF through this workflow, including structural variants.

References:
- `GATK-SV Package`: [v1.0.5](https://github.com/broadinstitute/gatk), which is used in the base image that is used in the [docker](dockerfiles/Dockerfile.AnnotateSVAnnotate).
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the GATK-SV featured workspace.


## External Allele Frequency (AF)
TODO

