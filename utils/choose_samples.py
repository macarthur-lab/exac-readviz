"""
Utility methods for choosing bam samples to display for a particular ExAC variant.
"""

import os
import re
import time


def choose_samples(het_or_hom, alt_allele_index, genotypes, sample_id_include_status):
    """Contains heuristics for choosing which samples to display for a given
    variant (chrom, pos, ref, alt) in it's het or hom-alt state.

    Args:
        het_or_hom: Either "het" or "hom" to indicate whether to choose het or hom-alt samples.
        alt_allele_index: 1-based index of the alt allele that we're choosing
            samples for (bewteen 1 and n = total number of alt alleles at this site)
        genotypes: a dictionary that maps each VCF sample_id to a 4-tuple: (gt_ref, gt_alt, GQ, DP)
            where gt_ref and gt_alt are integers between 0 and num alt alleles.
        sample_id_include_status: dictionary mapping each VCF sample_id to either
            True or False depending on the last column of the exac info table

    Return:
        Iterator over sample ids in order - with the ones that should
        displayed first.
    """
    assert het_or_hom in ["het", "hom"], "Unexpected het_or_hom arg: %s" % het_or_hom

    # filter all samples down to just samples that have the desired genotype and have include=YES
    relevant_samples = []  # a list of dicts
    for sample_id, (gt_ref, gt_alt, GQ, DP) in genotypes.items():
        if gt_ref is None and gt_alt is None:
            continue

        if het_or_hom == "het":
            if gt_ref == gt_alt:
                continue  # skip if homozygous
            if gt_ref != alt_allele_index and gt_alt != alt_allele_index:
                continue # skip if neither allele matches the specific alt allele we're looking for (eg. 1/3)
        elif het_or_hom == "hom":
            if gt_ref != gt_alt:
                continue  # skip unless homozygous
            if gt_alt != alt_allele_index:
                continue  # skip unless homozygous for the specific alt allele we're looking for (this matters for multiallelics)
        else:
            raise ValueError("Unexpected het_or_hom value: " + het_or_hom)

        if DP < 10 or GQ < 20:
            continue  # skip samples that don't pass _Adj thresholds since they are not counted in the ExAC browser het/hom counts.

        if not sample_id_include_status[sample_id]:
            continue  # skip samples where include status != "YES"

        relevant_samples.append( {"sample_id": sample_id, "GQ": GQ} )

    # return samples in order from highest to lowest GQ.
    for sample in sorted(relevant_samples, key=lambda s: s["GQ"], reverse=True):
        yield sample["sample_id"]  #, sample["GQ"]
