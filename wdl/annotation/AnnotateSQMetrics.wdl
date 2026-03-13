version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSQMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call CalculateSiteMetrics {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.site_quality",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs as MergeAnnotatedVcfs {
        input:
            vcfs = CalculateSiteMetrics.annotated_vcf,
            vcf_idxs = CalculateSiteMetrics.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.site_quality",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File site_quality_annotated_vcf = MergeAnnotatedVcfs.concat_vcf
        File site_quality_annotated_vcf_idx = MergeAnnotatedVcfs.concat_vcf_idx
    }
}

task CalculateSiteMetrics {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x INFO/AS_pab_max,INFO/AS_VarDP,INFO/AS_QUALapprox,INFO/AS_QD,INFO/inbreeding_coeff,INFO/HWE \
            -Oz -o stripped.vcf.gz \
            "~{vcf}"
        tabix -p vcf stripped.vcf.gz

        python3 <<CODE
import pysam
import scipy.stats as stats

vcf_in = pysam.VariantFile("stripped.vcf.gz")

vcf_in.header.info.add("AS_pab_max", "1", "Float", "Allele-specific max p-value for binomial test of allele balance (expected 0.5)")
vcf_in.header.info.add("AS_VarDP", "1", "Integer", "Allele-specific depth over variant genotypes")
vcf_in.header.info.add("AS_QUALapprox", "1", "Integer", "Allele-specific sum of PL[0] values; used to approximate the QUAL score")
vcf_in.header.info.add("AS_QD", "1", "Float", "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant")
vcf_in.header.info.add("inbreeding_coeff", "1", "Float", "Inbreeding coefficient from Hardy-Weinberg expectation")
vcf_in.header.info.add("HWE", "1", "Float", "Hardy-Weinberg equilibrium p-value")

vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'wz', header=vcf_in.header)

for rec in vcf_in:
    # Skip multiallelics
    if len(rec.alts) != 1:
        vcf_out.write(rec)
        continue

    vardp = 0
    qualapprox = 0
    pabs = []
    n_ref = n_het = n_alt = 0

    for sample in rec.samples.values():
        gt = sample.get('GT')
        if not gt or None in gt:
            continue

        # Count alt alleles to classify genotype (handles phased and unphased)
        n_alts = sum(gt)
        if n_alts == 0:
            n_ref += 1
        elif n_alts == 1:
            n_het += 1
        elif n_alts == 2:
            n_alt += 1

        # AS_VarDP: sum DP across variant carriers
        if n_alts > 0:
            dp = sample.get('DP')
            if dp is not None:
                vardp += dp
            pl = sample.get('PL')
            if pl and len(pl) > 0 and pl[0] is not None:
                qualapprox += pl[0]

        # AS_pab_max: binomial test on AD for heterozygotes
        if n_alts == 1 and 'AD' in sample:
            ad = sample['AD']
            if ad and len(ad) == 2 and ad[0] is not None and ad[1] is not None:
                total_ad = ad[0] + ad[1]
                if total_ad > 0:
                    pabs.append(stats.binomtest(ad[1], total_ad, 0.5).pvalue)

    # Set site-level metrics only when data is available
    if vardp > 0:
        rec.info['AS_VarDP'] = vardp
    if qualapprox > 0:
        rec.info['AS_QUALapprox'] = qualapprox
    if vardp > 0 and qualapprox > 0:
        rec.info['AS_QD'] = qualapprox / vardp
    if pabs:
        rec.info['AS_pab_max'] = max(pabs)

    n_tot = n_ref + n_het + n_alt
    if n_tot > 0:
        p = (2.0 * n_ref + n_het) / (2.0 * n_tot)
        q = 1.0 - p
        e_het = 2.0 * p * q * n_tot

        # Inbreeding coefficient: 1 - (observed hets / expected hets)
        if e_het > 0:
            rec.info['inbreeding_coeff'] = 1.0 - (n_het / e_het)

        # HWE p-value via chi-square goodness of fit
        e_ref = (p ** 2) * n_tot
        e_alt = (q ** 2) * n_tot
        chisq = sum(((o - e) ** 2) / e for o, e in zip([n_ref, n_het, n_alt], [e_ref, e_het, e_alt]) if e > 0)
        rec.info['HWE'] = stats.chi2.sf(chisq, 1)

    vcf_out.write(rec)

vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
