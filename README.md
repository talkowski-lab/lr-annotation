# gnomAD Long-Read
Methods and tools for callset analysis, QC and post-processing for the processing of gnomAD Long-Read data.


## Allele Frequency Annotations
We used the `AnnotateVcf` workflow from GATK-SV, albeit using a modified version from the `kj_gnomad_lr` branch - details regarding this modification can be seen in the following section.


## Functional Annotations
For SNVs, TODO.

For structural variants, we annotated them using a modified version of the `AnnotateVcf` workflow from the `kj_gnomad_lr` branch of GATK-SV. This version passes a new input `min_annotation_size=50` that indicates a minimum size threshold (based on the `INFO/SVLEN` field) to annotate variants for with `SVAnnotate`.


##  STR Annotations
TODO.


## Mobile Element Insertion (MEI) Annotations
TODO.


## Outstanding TODOs
- Update `AnnotateVcf` inputs to point to gnomAD V4 reference files as found in the featured workspace.
