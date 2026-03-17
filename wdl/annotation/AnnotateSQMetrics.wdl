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
            -x INFO/AS_VarDP,INFO/AS_QUALapprox,INFO/AS_QD,INFO/inbreeding_coeff,INFO/HWE,INFO/AS_pab_max \
            -Oz -o stripped.vcf.gz \
            "~{vcf}"
        tabix -p vcf stripped.vcf.gz

        python3 <<CODE
import pysam
import scipy.stats as stats

vcf_in = pysam.VariantFile("stripped.vcf.gz")

vcf_in.header.info.add("inbreeding_coeff", "A", "Float", "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium)")
vcf_in.header.info.add("AS_pab_max", "A", "Float", "Allele-specific maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5")
vcf_in.header.info.add("AS_QUALapprox", "A", "Integer", "Allele-specific sum of PL[0] values; used to approximate the QUAL score")
vcf_in.header.info.add("AS_QD", "A", "Float", "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant")
vcf_in.header.info.add("AS_VarDP", "A", "Integer", "Allele-specific depth over variant genotypes (does not include depth of reference samples)")
vcf_in.header.info.add("HWE", "1", "Float", "Hardy-Weinberg equilibrium p-value")

vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'wz', header=vcf_in.header)

for rec in vcf_in:
    if not rec.alts:
        vcf_out.write(rec)
        continue

    n_alt_alleles = len(rec.alts)
    as_vardp = [0] * n_alt_alleles
    as_qualapprox = [0] * n_alt_alleles
    as_qd_sum = [0.0] * n_alt_alleles
    as_qd_n = [0] * n_alt_alleles
    n_ref_per_allele = [0] * n_alt_alleles
    n_het_per_allele = [0] * n_alt_alleles
    n_alt_per_allele = [0] * n_alt_alleles
    n_ref_site = n_het_site = n_alt_site = 0
    as_pab_max = [None] * n_alt_alleles

    for sample in rec.samples.values():
        gt = sample.get('GT')
        if not gt or None in gt:
            continue

        n_non_ref = sum(1 for allele in gt if allele is not None and allele > 0)
        if n_non_ref == 0:
            n_ref_site += 1
        elif n_non_ref == 1:
            n_het_site += 1
        else:
            n_alt_site += 1

        pl = sample.get('PL')
        hom_ref_pl = None
        if pl and len(pl) > 0 and pl[0] is not None:
            hom_ref_pl = int(pl[0])

        dp = sample.get('DP')
        ad = sample.get('AD')
        ref_copies = sum(1 for allele in gt if allele == 0)

        for allele_idx in range(1, n_alt_alleles + 1):
            allele_i = allele_idx - 1
            n_copies = sum(1 for allele in gt if allele == allele_idx)

            if n_copies == 0:
                n_ref_per_allele[allele_i] += 1
            elif n_copies == 1:
                n_het_per_allele[allele_i] += 1
            else:
                n_alt_per_allele[allele_i] += 1

            if n_copies > 0:
                if dp is not None:
                    as_vardp[allele_i] += int(dp)
                if hom_ref_pl is not None:
                    as_qualapprox[allele_i] += hom_ref_pl
                    if ad and len(ad) > allele_idx and ad[allele_idx] is not None and ad[allele_idx] > 0:
                        as_qd_sum[allele_i] += hom_ref_pl / float(ad[allele_idx])
                        as_qd_n[allele_i] += 1

            if n_copies == 1 and ref_copies == 1 and ad is not None and len(ad) > allele_idx:
                ref_dp = ad[0]
                alt_dp = ad[allele_idx]
                if ref_dp is not None and alt_dp is not None:
                    total_dp = ref_dp + alt_dp
                    if total_dp > 0:
                        pval = round(float(stats.binomtest(ref_dp, total_dp, 0.5).pvalue), 6)
                        if as_pab_max[allele_i] is None or pval > as_pab_max[allele_i]:
                            as_pab_max[allele_i] = pval

    rec.info['AS_VarDP'] = tuple(as_vardp)
    rec.info['AS_QUALapprox'] = tuple(as_qualapprox)
    rec.info['AS_QD'] = tuple(
        round(as_qd_sum[i] / as_qd_n[i], 6) if as_qd_n[i] > 0 else 0.0
        for i in range(n_alt_alleles)
    )

    inbreeding_coeff = []
    for i in range(n_alt_alleles):
        n_ref_allele = n_ref_per_allele[i]
        n_het_allele = n_het_per_allele[i]
        n_alt_allele = n_alt_per_allele[i]
        n_tot_allele = n_ref_allele + n_het_allele + n_alt_allele

        if n_tot_allele == 0:
            inbreeding_coeff.append(0.0)
            continue

        p_allele = (2.0 * n_ref_allele + n_het_allele) / (2.0 * n_tot_allele)
        q_allele = 1.0 - p_allele
        e_het_allele = 2.0 * p_allele * q_allele * n_tot_allele

        if e_het_allele > 0:
            inbreeding_coeff.append(round(1.0 - (n_het_allele / e_het_allele), 6))
        else:
            inbreeding_coeff.append(0.0)

    rec.info['inbreeding_coeff'] = tuple(inbreeding_coeff)

    if any(v is not None for v in as_pab_max):
        rec.info['AS_pab_max'] = tuple(as_pab_max)

    n_tot = n_ref_site + n_het_site + n_alt_site
    if n_tot > 0:
        p = (2.0 * n_ref_site + n_het_site) / (2.0 * n_tot)
        q = 1.0 - p
        e_het = 2.0 * p * q * n_tot
        e_ref = (p ** 2) * n_tot
        e_alt = (q ** 2) * n_tot
        chisq = sum(((o - e) ** 2) / e for o, e in zip([n_ref_site, n_het_site, n_alt_site], [e_ref, e_het, e_alt]) if e > 0)
        rec.info['HWE'] = round(float(stats.chi2.sf(chisq, 1)), 6)

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
