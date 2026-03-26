version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillSVFormatFields {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[File] sample_sv_stats
        String prefix

        String utils_docker

        File? swap_samples

        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_compute_counts
        RuntimeAttr? runtime_attr_aggregate_counts
    }

    if (defined(swap_samples)) {
        call Helpers.SwapSampleIds {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                sample_swap_list = select_first([swap_samples]),
                prefix = "~{prefix}.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples
        }
    }

    File final_cohort_vcf = select_first([SwapSampleIds.swapped_vcf, cohort_vcf])
    File final_cohort_vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, cohort_vcf_idx])

    scatter (i in range(length(sample_ids))) {
        call ComputePerSampleCallerCounts {
            input:
                sample_id = sample_ids[i],
                sv_stats = sample_sv_stats[i],
                prefix = "~{prefix}.~{sample_ids[i]}.caller_counts",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_compute_counts
        }
    }

    call AggregateCallerCounts {
        input:
            per_sample_counts = ComputePerSampleCallerCounts.counts_tsv,
            prefix = "~{prefix}.caller_counts",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_aggregate_counts
    }

    output {
        File out_cohort_vcf = final_cohort_vcf
        File out_cohort_vcf_idx = final_cohort_vcf_idx
        File caller_counts_tsv = AggregateCallerCounts.counts_tsv
    }
}

task ComputePerSampleCallerCounts {
    input {
        String sample_id
        File sv_stats
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip

def get_len_bucket(svlen):
    l = abs(int(svlen))
    if l < 50:
        return "<50"
    elif l < 500:
        return "50-500"
    elif l < 5000:
        return "500-5000"
    else:
        return ">5000"

counts = {}

with gzip.open("~{sv_stats}", 'rt') as f:
    header = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if header is None:
            header = line.lstrip('#').split('\t')
            continue
        fields = dict(zip(header, line.split('\t')))
        svtype = fields['SVTYPE']
        if svtype not in ('INS', 'DEL'):
            continue
        len_bucket = get_len_bucket(fields['SVLEN'])
        for caller in fields['SUPP'].split(','):
            caller = caller.strip()
            key = (svtype, len_bucket, caller)
            counts[key] = counts.get(key, 0) + 1

with open("~{prefix}.tsv", 'w') as out:
    out.write("TYPE\tLEN\tCALLER\tCOUNT\n")
    for (svtype, len_bucket, caller), count in sorted(counts.items()):
        out.write(f"{svtype}\t{len_bucket}\t{caller}\t{count}\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(sv_stats, "GB")) + 5,
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

task AggregateCallerCounts {
    input {
        Array[File] per_sample_counts
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv
from collections import defaultdict

input_files = "~{sep=',' per_sample_counts}".split(',')

agg = defaultdict(int)
for f in input_files:
    with open(f, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = (row['TYPE'], row['LEN'], row['CALLER'])
            agg[key] += int(row['COUNT'])

callers = sorted(set(k[2] for k in agg.keys()))
len_order = ["<50", "50-500", "500-5000", ">5000"]
type_order = ["INS", "DEL"]

with open("~{prefix}.tsv", 'w') as out:
    header = ["TYPE", "LEN"] + [f"COUNT_{c}" for c in callers]
    out.write("\t".join(header) + "\n")
    for svtype in type_order:
        for len_bucket in len_order:
            row = [svtype, len_bucket] + [str(agg.get((svtype, len_bucket, c), 0)) for c in callers]
            out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(per_sample_counts, "GB")) + 5,
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
