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
                prefix = prefix + "_shard_" + shard_idx
        }

        call Aln.AlignGeneric as alignIns{
            input:
                input_fa = InsToFaWithFlanking.ins_fa,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = prefix + "_shard_" + shard_idx,
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

     call Helpers.MergeVcfs { input: 
        vcfs=PALMER.vcf, 
        vcf_idxs = PALMER.vcf_tbi,
        prefix=prefix+"_merged",
        docker=utils_docker}

    output {
        File PALMER_vcf = MergeVcfs.merged_vcf
        File PALMER_vcf_tbi = MergeVcfs.merged_vcf_idx
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
    left_flank = ref_fa.fetch(chrom, max(0,pos-10000), pos)
    right_flank = ref_fa.fetch(chrom, pos, min(pos+10000, chrom_length))
    outfile.write('>%s_%d_%s_%s_%d\n%s\n' % (chrom, pos, ref, ID, counter, left_flank+ins_seq[1:]+right_flank))
infile.close()
outfile.close()
EOF
    >>>

    output {
        File ins_fa = '~{prefix}.fa'
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
        bcftools view -h ~{prefix}_integrated.vcf > header.txt
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">' >> header.txt
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