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
        String MEI_type
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

        call PALMER {
            input:
                bam = alignIns.bamOut,
                bai = alignIns.baiOut,
                ref_fa = ref_fa,
                prefix = prefix,
                mode = mode,
                MEI_type = MEI_type,
        }
     }

    call Helpers.MergeVcfs { 
        input: 
            vcfs=PALMER.vcf, 
            vcf_idxs = PALMER.vcf_tbi,
            prefix=prefix+"_merged",
            docker=utils_docker,
            extra_args="--force-samples"
    }

    call Helpers.ConcatTsvs as ConcatTSDReads {
        input:
            tsvs = PALMER.TSD_reads,
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
            ins_variant_keys = ConcatInsVariantKeys.concatenated_tsv,
            prefix = prefix
    }

    output {
        File PALMER_vcf = MergeVcfs.merged_vcf
        File PALMER_vcf_tbi = MergeVcfs.merged_vcf_idx
        File TSD_reads = ConcatTSDReads.concatenated_tsv
        File annotated_vcf = AnnotateOriginalVcfWithPALMER.annotated_vcf
        File annotated_vcf_tbi = AnnotateOriginalVcfWithPALMER.annotated_vcf_tbi
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
ref_fa = FastaFile("~{ref_fa}")
infile = open("~{prefix}.tmp.txt", "r")
outfile = open("~{prefix}.fa", "w")
mapfile = open("~{prefix}.variant_keys.tsv", "w")
flanking_bp = ~{flanking_bp}

prev_chrom = ''
counter=0 #counter is to guarantee that no seqs have identical names in the output fasta
chrom_length=0

for line in infile:
    chrom, pos, ID, ref, ins_seq = line.strip().split('\t')
    pos = int(pos)
    if prev_chrom!=chrom:
        counter=0
        prev_chrom=chrom
        chrom_length = ref_fa.get_reference_length(chrom)
    else:
        counter+=1
    left_flank = ref_fa.fetch(chrom, max(0, pos-flanking_bp), pos)
    right_flank = ref_fa.fetch(chrom, pos, min(pos+flanking_bp, chrom_length))
    fasta_name = '%s_%d_%s_%s_%d' % (chrom, pos, ref, ID, counter)
    outfile.write('>%s\n%s\n' % (fasta_name, left_flank+ins_seq[1:]+right_flank))
    mapfile.write('%s\t%s\t%d\t%s\t%s\t%s\n' % (fasta_name, chrom, pos, ID, ref, (ref + ins_seq[1:]).upper()))
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

        #fix output VCF header (missing SVTYPE definition)
        bcftools view -h ~{prefix}_integrated.vcf |grep -v '^#CHROM' > header.txt
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">' >> header.txt
        bcftools view -h ~{prefix}_integrated.vcf|tail -n1 >> header.txt
        bcftools reheader -h header.txt ~{prefix}_integrated.vcf > ~{prefix}_integrated_fixed.vcf
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
        File ins_variant_keys
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
from collections import OrderedDict

import pysam


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

def normalize_scalar(value):
    if isinstance(value, tuple):
        if len(value) == 0:
            return None
        if len(value) == 1:
            return value[0]
        return ",".join(str(v) for v in value)
    return value


def is_missing(value):
    return value is None or value == "." or value == ""


def get_subtype_label(record):
    subtype = normalize_scalar(record.info.get("SUBTYPE"))
    if not is_missing(subtype) and subtype != "NA":
        return str(subtype)
    return "UNKNOWN"


header_to_key = {}
with open("~{ins_variant_keys}", "r", encoding="utf-8") as map_handle:
    for raw_line in map_handle:
        raw_line = raw_line.rstrip("\n")
        if not raw_line:
            continue
        fasta_name, chrom, pos_str, variant_id, ref, alt = raw_line.split("\t")
        pos = int(pos_str)
        header_to_key[fasta_name] = (chrom, pos, variant_id, ref.upper(), alt.upper())

palmer_annotation_candidates = {}
palmer_vcf = pysam.VariantFile("~{palmer_vcf}")
for record in palmer_vcf:
    if "PASS" not in set(record.filter.keys()):
        continue
    source_name = normalize_scalar(record.info.get("TSD_READ"))
    if is_missing(source_name) or source_name == "NA":
        continue
    key = header_to_key.get(source_name)
    if key is None:
        continue

    annotation = {}
    for field, (source_field, _field_type, _description) in TRANSFER_FIELDS.items():
        if field == "PALMER_MEI_LEN":
            continue
        value = normalize_scalar(record.info.get(source_field))
        if is_missing(value) or value == "NA":
            continue
        annotation[field] = value

    ins_seq = normalize_scalar(record.info.get("INS_SEQ"))
    if not is_missing(ins_seq) and ins_seq != "NA":
        annotation["PALMER_MEI_LEN"] = len(str(ins_seq))

    candidate_info = palmer_annotation_candidates.setdefault(
        key,
        {"annotations": [], "subtypes": []},
    )
    if annotation not in candidate_info["annotations"]:
        candidate_info["annotations"].append(annotation)
        candidate_info["subtypes"].append(get_subtype_label(record))

palmer_vcf.close()

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

original_vcf = pysam.VariantFile("~{original_vcf}")
header = original_vcf.header.copy()
for field, (_source_field, field_type, description) in TRANSFER_FIELDS.items():
    if field not in header.info:
        header.add_line(
            f'##INFO=<ID={field},Number=1,Type={field_type},Description="{description}">'
        )
multi_field, (multi_field_type, multi_field_description) = MULTI_MEI_FIELD
if multi_field not in header.info:
    header.add_line(
        f'##INFO=<ID={multi_field},Number=.,Type={multi_field_type},Description="{multi_field_description}">'
    )

annotated_count = 0
multi_mei_annotated_count = 0
with pysam.VariantFile("~{prefix}.palmer_annotated.vcf", "w", header=header) as out_vcf:
    for record in original_vcf:
        matches = []
        record_id = record.id if record.id is not None else "."
        if record.alts:
            for alt in record.alts:
                key = (record.contig, record.pos, record_id, record.ref.upper(), alt.upper())
                if key in palmer_annotations:
                    matches.append(key)

        if len(matches) > 1:
            raise RuntimeError(
                "Multiple PALMER matches found for original record "
                + f"{record.contig}:{record.pos}:{record_id}"
            )

        if len(matches) == 1:
            annotation = palmer_annotations[matches[0]]
            for field, value in annotation.items():
                record.info[field] = value
            annotated_count += 1
            if multi_field in annotation:
                multi_mei_annotated_count += 1

        out_vcf.write(record)

original_vcf.close()

with open("~{prefix}.palmer_annotation_summary.txt", "w", encoding="utf-8") as summary:
    summary.write(f"annotated_records\t{annotated_count}\n")
    summary.write(f"annotated_multi_mei_records\t{multi_mei_annotated_count}\n")
    summary.write(f"palmer_multi_mei_source_keys\t{multi_mei_count}\n")
    summary.write(f"palmer_records_with_source_key\t{len(palmer_annotations)}\n")
EOF

        bcftools sort -Oz -o ~{prefix}.palmer_annotated.vcf.gz ~{prefix}.palmer_annotated.vcf
        tabix -p vcf ~{prefix}.palmer_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.palmer_annotated.vcf.gz"
        File annotated_vcf_tbi = "~{prefix}.palmer_annotated.vcf.gz.tbi"
        File annotation_summary = "~{prefix}.palmer_annotation_summary.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(original_vcf, "GB") + size(palmer_vcf, "GB") + size(ins_variant_keys, "GB")) + 20,
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