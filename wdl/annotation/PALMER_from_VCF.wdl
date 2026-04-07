version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../tools/MinimapAlignment.wdl" as Aln

workflow CallPALMER {
    input {
        File vcf
        File vcf_tbi
        File ref_fa
        File ref_fai
        String prefix
        String mode
        Array[String] MEI_types = ["ALU","LINE","SVA"]
        String utils_docker="us-central1-docker.pkg.dev/talkowski-training/kj-development/utils:kj_V9"

        Int records_per_shard = 50000
        Int flanking_bp = 10000
    }

    call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_tbi,
                records_per_shard = records_per_shard,
                prefix = prefix,
                docker = utils_docker
    }

    scatter (shard_idx in range(length(ShardVcfByRecords.shards))) {
        call InsToFaWithFlanking {
            input:
                vcf = ShardVcfByRecords.shards[shard_idx],
                vcf_tbi = ShardVcfByRecords.shard_idxs[shard_idx],
                ref_fa = ref_fa,
                flanking_bp = flanking_bp,
                prefix = "~{prefix}_shard_~{shard_idx}"
        }

        call Aln.AlignGeneric as alignIns{
            input:
                input_fa = InsToFaWithFlanking.ins_fa,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}_shard_~{shard_idx}",
                docker = "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28",
                flags = "-ayYL --MD -x asm5"
        }

        scatter (mei_type in MEI_types) {
            call PALMER {
                input:
                    bam = alignIns.bamOut,
                    bai = alignIns.baiOut,
                    ref_fa = ref_fa,
                    prefix = "~{prefix}_shard_~{shard_idx}_~{mei_type}",
                    mode = mode,
                    MEI_type = mei_type,
            }
        }
    }

    Array[File] palmer_vcfs = flatten(PALMER.vcf)
    Array[File] palmer_vcf_idxs = flatten(PALMER.vcf_tbi)
    Array[File] palmer_tsd_reads = flatten(PALMER.TSD_reads)

    call Helpers.MergeVcfs { 
        input: 
            vcfs=palmer_vcfs, 
            vcf_idxs = palmer_vcf_idxs,
            prefix=prefix+"_merged",
            docker=utils_docker,
            extra_args="--force-samples"
    }

    call Helpers.ConcatTsvs as ConcatTSDReads {
        input:
            tsvs = palmer_tsd_reads,
            prefix = prefix + "_TSD_reads",
            docker = utils_docker,
            preserve_header = true,
            skip_sort = true
    }

    call Helpers.ConcatTsvs as ConcatInsVariantKeys {
        input:
            tsvs = InsToFaWithFlanking.variant_key_tsv,
            prefix = prefix + "_ins_variant_keys",
            docker = utils_docker
    }

    call AnnotateOriginalVcfWithPALMER {
        input:
            original_vcf = vcf,
            original_vcf_idx = vcf_tbi,
            palmer_vcf = MergeVcfs.merged_vcf,
            palmer_vcf_idx = MergeVcfs.merged_vcf_idx,
            palmer_tsd_reads = ConcatTSDReads.concatenated_tsv,
            ins_variant_keys = ConcatInsVariantKeys.concatenated_tsv,
            prefix = prefix
    }

    output {
        File PALMER_vcf = MergeVcfs.merged_vcf
        File PALMER_vcf_tbi = MergeVcfs.merged_vcf_idx
        File TSD_reads = ConcatTSDReads.concatenated_tsv
        File annotated_vcf = AnnotateOriginalVcfWithPALMER.annotated_vcf
        File annotated_vcf_tbi = AnnotateOriginalVcfWithPALMER.annotated_vcf_tbi
        File annotation_header = AnnotateOriginalVcfWithPALMER.annotation_header
        File annotation_tsv_gz = AnnotateOriginalVcfWithPALMER.annotation_tsv_gz
        File annotation_tsv_gz_tbi = AnnotateOriginalVcfWithPALMER.annotation_tsv_gz_tbi
    }
}

task InsToFaWithFlanking {
    input {
        File vcf
        File vcf_tbi
        Int flanking_bp=10000
        File ref_fa
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -i 'SVTYPE=="INS"' \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{vcf} > ~{prefix}.tmp.txt
        
        python3 <<EOF
import sys
from pysam import FastaFile

# Read the reference so each insertion can be exported with reference flanks.
ref_fa = FastaFile("~{ref_fa}")
infile = open("~{prefix}.tmp.txt", "r")
outfile = open("~{prefix}.fa", "w")
# Keep a lookup from the generated FASTA header back to the original variant
# identity used later when PALMER calls are mapped onto the input VCF.
mapfile = open("~{prefix}.variant_keys.tsv", "w")
flanking_bp = ~{flanking_bp}

prev_chrom = ''
# counter is included in the FASTA ID to avoid potential duplicate IDs
counter=0
chrom_length=0

for line in infile:
    chrom, pos, ID, ref, ins_seq = line.strip().split('\t')
    pos = int(pos)
    # Reset the per-chromosome counter when moving to a new contig
    # Save the contig length so right flanks do not extend past the reference.
    if prev_chrom!=chrom:
        counter=0
        prev_chrom=chrom
        chrom_length = ref_fa.get_reference_length(chrom)
    else:
        counter+=1

    # Extract reference sequence flanking the insertion breakpoint
    left_flank = ref_fa.fetch(chrom, max(0, pos-flanking_bp), pos)
    right_flank = ref_fa.fetch(chrom, pos, min(pos+flanking_bp, chrom_length))

    # Encode locus + allele context into the FASTA name so PALMER output can
    # later be traced back to the original input record.
    fasta_name = '%s_%d_%s_%s_%d' % (chrom, pos, ref, ID, counter)

    # Build the sequence given to PALMER: left flank + inserted bases + right
    # flank. ins_seq[1:] drops the anchor base because VCF INS alleles are
    # represented as REF + inserted sequence.
    outfile.write('>%s\n%s\n' % (fasta_name, left_flank+ins_seq[1:]+right_flank))

    # Write a mapping table with:
    #   generated FASTA name, chrom, pos, original ID, REF, ALT
    mapfile.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (fasta_name, chrom, pos, ID, ref, ins_seq))
infile.close()
outfile.close()
mapfile.close()
EOF
    >>>

    output {
        File ins_fa = '~{prefix}.fa'
        File variant_key_tsv = '~{prefix}.variant_keys.tsv'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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
        docker: "quay.io/ymostovoy/lr-process-mendelian:1.0"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PALMER {
    input {
        File bam
        File bai
        File ref_fa
        String prefix
        String mode
        String MEI_type

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "aligned bam, can contain long reads or assembled contigs"
        bai:              "index accompanying the BAM"
        ref_fa:           "fasta of reference used"
        prefix:           "output file prefix"
        mode:             "either raw or asm"
        MEI_type:         "type of MEI to call, options are LINE, ALU, SVA, HERVK"
    }

    Int disk_size = ceil(size(bam, "GB") + size(ref_fa, "GB")) * 6

    command <<<
        set -x

        mv ~{bam} ./
        mv ~{bai} ./
        
        bamfile=$(basename ~{bam})

        /PALMER/PALMER --input ${bamfile} \
                --ref_fa ~{ref_fa} \
                --ref_ver GRCh38 \
                --type ~{MEI_type} \
                --mode ~{mode} \
                --output "~{prefix}" \
                --chr ALL \
                --workdir ./

        python3 <<EOF
import csv

mei_type = "~{MEI_type}"
# Prefix all PALMER cluster IDs with the MEI type so IDs remain unique
# after the workflow scatters over multiple MEI classes and merges outputs.
cluster_prefix = f"{mei_type}_"

def prefix_text_table(path):
    # Rewrite the first column of PALMER's tabular outputs so the cluster IDs
    # match the MEI-prefixed IDs that will also be written into the VCF.
    rows = []
    with open(path, "r", encoding="utf-8") as src:
        reader = csv.reader(src, delimiter="\t")
        for row in reader:
            if not row:
                rows.append(row)
                continue
            # Keep the header unchanged; prefix only data rows.
            if row[0] != "cluster_id":
                row[0] = cluster_prefix + row[0]
            rows.append(row)

    with open(path, "w", encoding="utf-8", newline="") as dst:
        writer = csv.writer(dst, delimiter="\t")
        writer.writerows(rows)

# Update PALMER's VCF IDs to use the same prefix as the tabular outputs.
# This prevents collisions between e.g. LINE/ALU/SVA cluster IDs after merge.
vcf_in = "~{prefix}_integrated.vcf"
vcf_out = "~{prefix}_integrated_namespaced.vcf"
with open(vcf_in, "r", encoding="utf-8") as src, open(vcf_out, "w", encoding="utf-8") as dst:
    for line in src:
        if line.startswith("#"):
            dst.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        # Column 3 is the VCF record ID. Leave missing IDs alone.
        if len(fields) > 2 and fields[2] not in (".", ""):
            fields[2] = cluster_prefix + fields[2]
        dst.write("\t".join(fields) + "\n")

# Apply the same ID rewrite to PALMER's auxiliary text reports so downstream
# joins against cluster IDs stay consistent with the rewritten VCF.
prefix_text_table("~{prefix}_TSD_reads_output.txt")
prefix_text_table("~{prefix}_calls.txt")
prefix_text_table("~{prefix}_all_reads_output.txt")
EOF

        #fix output VCF header (missing SVTYPE definition)
        bcftools view -h ~{prefix}_integrated_namespaced.vcf |grep -v '^#CHROM' > header.txt
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">' >> header.txt
        bcftools view -h ~{prefix}_integrated_namespaced.vcf|tail -n1 >> header.txt
        bcftools reheader -h header.txt ~{prefix}_integrated_namespaced.vcf > ~{prefix}_integrated_fixed.vcf
        bcftools sort -Oz -o ~{prefix}_integrated.vcf.gz ~{prefix}_integrated_fixed.vcf
        tabix -p vcf ~{prefix}_integrated.vcf.gz

    >>>

    output {
        File vcf = "~{prefix}_integrated.vcf.gz"
        File vcf_tbi = "~{prefix}_integrated.vcf.gz.tbi"
        File TSD_reads = "~{prefix}_TSD_reads_output.txt"
        File all_reads = "~{prefix}_all_reads_output.txt"
        File calls = "~{prefix}_calls.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-palmer:2.3.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task AnnotateOriginalVcfWithPALMER {
    input {
        File original_vcf
        File original_vcf_idx
        File palmer_vcf
        File palmer_vcf_idx
        File palmer_tsd_reads
        File ins_variant_keys
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
from collections import OrderedDict

import pysam

# Mapping of PALMER INFO fields to copy onto the original VCF, along with the header type/description
TRANSFER_FIELDS = OrderedDict([
    ("PALMER_MEI_TYPE", ("SUBTYPE", "String", "PALMER MEI subtype")),
    ("ORIENTATION", ("ORIENTATION", "String", "Insertion orientation")),
    ("POLYA", ("POLYA", "Integer", "polyA tail size")),
    ("TSD5", ("TSD5", "Integer", "5' TSD size")),
    ("TSD3", ("TSD3", "Integer", "3' TSD size")),
    ("TRANSD", ("TRANSD", "Integer", "Predicted 3' transduction size")),
    ("INV5", ("INV5", "String", "Has 5' inverted sequence")),
    ("INV5_START", ("INV5_START", "Integer", "Start position of 5' inverted sequence")),
    ("INV5_END", ("INV5_END", "Integer", "End position of 5' inverted sequence")),
    ("START_INVAR", ("START_INVAR", "Integer", "Average start position within insertion sequence")),
    ("END_INVAR", ("END_INVAR", "Integer", "Average end position within insertion sequence")),
    ("TSD5_SEQ", ("TSD5_SEQ", "String", "5' TSD sequence (NA if unavailable)")),
    ("TSD3_SEQ", ("TSD3_SEQ", "String", "3' TSD sequence (NA if unavailable)")),
    ("TRANSD_READ", ("TRANSD_READ", "String", "Predicted transduction sequence (NA if unavailable)")),
    ("JUNC_26MER", ("JUNC_26MER", "String", "Unique 26-mer at 5' junction (NA if unavailable)")),
    ("PALMER_MEI_LEN", ("INS_SEQ", "Integer", "Length of PALMER INS_SEQ annotation")),
])
MULTI_MEI_FIELD = (
    "PALMER_multi_MEI",
    ("String", "Multiple PALMER calls detected for this insertion; other PALMER INFO annotations were not transferred"),
)
ANNOTATION_FIELDS = list(TRANSFER_FIELDS.keys()) + [MULTI_MEI_FIELD[0]]

def normalize_scalar(value):
    # Pysam may return tuples for INFO fields; collapse them to a comparable
    # scalar representation so values can be stored and deduplicated easily.
    if isinstance(value, tuple):
        if len(value) == 0:
            return None
        if len(value) == 1:
            return value[0]
        return ",".join(str(v) for v in value)
    return value


def is_missing(value):
    # Treat common VCF missing-value conventions as absent.
    return value is None or value == "." or value == ""


def get_subtype_label(record):
    # Use SUBTYPE when available; otherwise retain a placeholder for the
    # multi-MEI fallback annotation.
    subtype = normalize_scalar(record.info.get("SUBTYPE"))
    if not is_missing(subtype) and subtype != "NA":
        return str(subtype)
    return "UNKNOWN"


# Load the lookup written by `InsToFaWithFlanking`, which maps each generated
# FASTA header back to the exact original VCF allele key.
fastaID_to_key = {}
with open("~{ins_variant_keys}", "r", encoding="utf-8") as map_handle:
    for raw_line in map_handle:
        raw_line = raw_line.rstrip("\n")
        if not raw_line:
            continue
        fasta_name, chrom, pos_str, variant_id, ref, alt = raw_line.split("\t")
        pos = int(pos_str)
        fastaID_to_key[fasta_name] = (chrom, pos, variant_id, ref, alt)

# PALMER's TSD-read output links a cluster ID to one or more source FASTA names.
# Capture those links so PALMER VCF calls can be traced back to the input VCF.
cluster_to_source_names = {}
with open("~{palmer_tsd_reads}", "r", encoding="utf-8") as tsd_handle:
    for raw_line in tsd_handle:
        raw_line = raw_line.rstrip("\n")
        if not raw_line or raw_line.startswith("cluster_id\t"):
            continue
        fields = raw_line.split("\t")
        if len(fields) < 2:
            continue
        cluster_id = fields[0]
        source_name = fields[1]
        cluster_to_source_names.setdefault(cluster_id, set()).add(source_name)

# Collect MEI-annotated VCF insertion keys and their corresponding PALMER annotation fields
palmer_annotation_candidates = {}
palmer_vcf = pysam.VariantFile("~{palmer_vcf}")
for record in palmer_vcf:
    # Only transfer annotations from calls that passed PALMER filtering.
    if "PASS" not in set(record.filter.keys()):
        continue

    # Source FASTA names can come from either the cluster-to-read table or the
    # representative `TSD_READ` INFO field embedded in the PALMER VCF record
    representative_source_name = normalize_scalar(record.info.get("TSD_READ"))
    source_names = set()
    if not is_missing(record.id):
        source_names.update(cluster_to_source_names.get(record.id, set()))
    if not is_missing(representative_source_name) and representative_source_name != "NA":
        source_names.add(representative_source_name)
    if not source_names:
        continue

    # Extract the PALMER annotations that should be added back to the original VCF record
    annotations = {}
    for field, (source_field, _field_type, _description) in TRANSFER_FIELDS.items():
        if field == "PALMER_MEI_LEN":
            continue
        value = normalize_scalar(record.info.get(source_field))
        if is_missing(value) or value == "NA":
            continue
        annotations[field] = value

    # Record the inserted sequence length
    ins_seq = normalize_scalar(record.info.get("INS_SEQ"))
    if not is_missing(ins_seq) and ins_seq != "NA":
        annotations["PALMER_MEI_LEN"] = len(str(ins_seq))

    for source_name in source_names:
        # Convert each PALMER source name back into the original VCF key
        key = fastaID_to_key.get(source_name)
        if key is None:
            continue
        candidate_info = palmer_annotation_candidates.setdefault(
            key, {"annotations": [], "subtypes": []},)

        # Avoid exact duplicate annotations but otherwise transfer annos for each original VCF key
        if annotations not in candidate_info["annotations"]:
            candidate_info["annotations"].append(annotations)
            candidate_info["subtypes"].append(get_subtype_label(record))

palmer_vcf.close()

# Resolve the candidate PALMER annotations for each original insertion. If more
# than one distinct PALMER call maps back to the same input record, emit only a
# summary field listing the competing subtypes
palmer_annotations = {}
multi_mei_count = 0
for key, candidate_info in palmer_annotation_candidates.items():
    annotations = candidate_info["annotations"]
    if len(annotations) <= 1:
        if annotations:
            palmer_annotations[key] = annotations[0]
        continue

    palmer_annotations[key] = {
        MULTI_MEI_FIELD[0]: tuple(candidate_info["subtypes"])
    }
    multi_mei_count += 1

# Emit new INFO header lines required by `bcftools annotate`.
multi_field, (multi_field_type, multi_field_description) = MULTI_MEI_FIELD
original_vcf = pysam.VariantFile("~{original_vcf}")
with open("~{prefix}.palmer_annotation_header.txt", "w", encoding="utf-8") as header_handle:
    for field, (_source_field, field_type, description) in TRANSFER_FIELDS.items():
        if field not in original_vcf.header.info:
            header_handle.write(
                f'##INFO=<ID={field},Number=1,Type={field_type},Description="{description}">\n'
            )
    if multi_field not in original_vcf.header.info:
        header_handle.write(
            f'##INFO=<ID={multi_field},Number=.,Type={multi_field_type},Description="{multi_field_description}">\n'
        )
original_vcf.close()

# Write annotations to a TSV to use with bcftools annotate
with open("~{prefix}.palmer_annotations.unsorted.tsv", "w", encoding="utf-8") as annotation_handle:
    for key in sorted(palmer_annotations):
        chrom, pos, variant_id, ref, alt = key
        annotation = palmer_annotations[key]
        row = [chrom, str(pos), variant_id, ref, alt]
        for field in ANNOTATION_FIELDS:
            value = annotation.get(field, ".")
            if isinstance(value, tuple):
                value = ",".join(str(v) for v in value)
            elif value is None or value == "":
                value = "."
            row.append(str(value))
        annotation_handle.write("\t".join(row) + "\n")

annotated_count = len(palmer_annotations)
multi_mei_annotated_count = sum(
    1 for annotation in palmer_annotations.values() if multi_field in annotation
)

# Emit a short summary to stdout to help verify how many records were annotated
# and how many required the multi-MEI fallback.
print(f"annotated_records\t{annotated_count}")
print(f"annotated_multi_mei_records\t{multi_mei_annotated_count}")
print(f"palmer_multi_mei_source_keys\t{multi_mei_count}")
print(f"palmer_records_with_source_key\t{len(palmer_annotations)}")
EOF

        LC_ALL=C sort -t $'\t' -k1,1 -k2,2n ~{prefix}.palmer_annotations.unsorted.tsv | bgzip -c > ~{prefix}.palmer_annotations.tsv.gz
        tabix -s 1 -b 2 -e 2 ~{prefix}.palmer_annotations.tsv.gz

        if [ -s ~{prefix}.palmer_annotations.tsv.gz ]; then
            if [ -s ~{prefix}.palmer_annotation_header.txt ]; then
                bcftools annotate \
                    -a ~{prefix}.palmer_annotations.tsv.gz \
                    -h ~{prefix}.palmer_annotation_header.txt \
                    -c CHROM,POS,~ID,REF,ALT,PALMER_MEI_TYPE,ORIENTATION,POLYA,TSD5,TSD3,TRANSD,INV5,INV5_START,INV5_END,START_INVAR,END_INVAR,TSD5_SEQ,TSD3_SEQ,TRANSD_READ,JUNC_26MER,PALMER_MEI_LEN,PALMER_multi_MEI \
                    -Oz -o ~{prefix}.palmer_annotated.vcf.gz \
                    ~{original_vcf}
            else # this is if the original VCF already had all necessary INFO headers, so no header file is needed for annotation
                bcftools annotate \
                    -a ~{prefix}.palmer_annotations.tsv.gz \
                    -c CHROM,POS,~ID,REF,ALT,PALMER_MEI_TYPE,ORIENTATION,POLYA,TSD5,TSD3,TRANSD,INV5,INV5_START,INV5_END,START_INVAR,END_INVAR,TSD5_SEQ,TSD3_SEQ,TRANSD_READ,JUNC_26MER,PALMER_MEI_LEN,PALMER_multi_MEI \
                    -Oz -o ~{prefix}.palmer_annotated.vcf.gz \
                    ~{original_vcf}
            fi
        else
            echo "No PALMER annotations to add; copying original VCF"
            cp ~{original_vcf} ~{prefix}.palmer_annotated.vcf.gz
        fi

        tabix -p vcf ~{prefix}.palmer_annotated.vcf.gz

        touch ~{prefix}.palmer_annotations.tsv.gz
        touch ~{prefix}.palmer_annotations.tsv.gz.tbi
        touch ~{prefix}.palmer_annotation_header.txt
    >>>

    output {
        File annotated_vcf = "~{prefix}.palmer_annotated.vcf.gz"
        File annotated_vcf_tbi = "~{prefix}.palmer_annotated.vcf.gz.tbi"
        File annotation_header = "~{prefix}.palmer_annotation_header.txt"
        File annotation_tsv_gz = "~{prefix}.palmer_annotations.tsv.gz"
        File annotation_tsv_gz_tbi = "~{prefix}.palmer_annotations.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(original_vcf, "GB") + size(palmer_vcf, "GB") + size(palmer_tsd_reads, "GB") + size(ins_variant_keys, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "quay.io/ymostovoy/lr-process-mendelian:1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: select_first([runtime_attr.docker, default_attr.docker])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}