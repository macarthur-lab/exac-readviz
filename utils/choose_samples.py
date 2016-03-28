"""
Utility methods for choosing bam samples to display for a particular ExAC variant.
"""

import os
import re
import time


def _is_hemizygous_segment(chrom, pos):
    return (chrom == 'X' and not ((60001 <= pos <= 2699520) or (154931044 <= pos <= 155260560))) or (chrom == 'Y')


def best_for_readviz_sample_id_iter(chrom, pos, het_or_hom_or_hemi, alt_allele_index, genotypes, sample_id_include_status, sample_id_sex):
    """Implements heuristics for choosing which samples are best to display for
    a given variant (chrom, pos, ref, alt) in it's het or hom-alt state.

    Args:
        chrom: chromosome string "chr1", "chrX", etc.
        pos:  variant position as an integer (minified)
        het_or_hom_or_hemi: Either "het" or "hom" or "hemi" to indicate whether to choose het, hom-alt, or hemi samples.
        alt_allele_index: 1-based index of the alt allele that we're choosing
            samples for (bewteen 1 and n = total number of alt alleles at this site)
        genotypes: a dictionary that maps each VCF sample_id to a 4-tuple: (gt_ref, gt_alt, GQ, DP)
            where gt_ref and gt_alt are integers between 0 and num alt alleles.
        sample_id_include_status: dictionary mapping each VCF sample_id to either
            True or False depending on the last column of the exac info table
        sample_id_sex: dictionary mapping each VCF sample_id to "m" or "f"

    Return:
        Iterator over sample ids in order - from those that should be first to
        be displayed to those that should be last.
    """
    assert het_or_hom_or_hemi in ["het", "hom", "hemi"], "Unexpected het_or_hom_or_hemi arg: %s" % het_or_hom_or_hemi

    # filter all samples down to just samples that have the desired genotype and have include=YES
    relevant_samples = []  # a list of dicts
    for sample_id, (gt_ref, gt_alt, GQ, DP) in genotypes.items():
        if gt_ref is None and gt_alt is None:
            continue

        if het_or_hom_or_hemi == "het":
            if gt_ref == gt_alt:
                continue  # skip if homozygous
            if gt_ref != alt_allele_index and gt_alt != alt_allele_index:
                continue # skip if neither allele matches the specific alt allele we're looking for (eg. 1/3)
        elif het_or_hom_or_hemi in ("hom", "hemi"):
            if gt_ref != gt_alt:
                continue  # skip unless homozygous
            if gt_alt != alt_allele_index:
                continue  # skip unless homozygous for the specific alt allele we're looking for (this matters for multiallelics)

            # handle hom vs. hemi or chromosomes X, Y
            if chrom not in ('X', 'Y'):
                if het_or_hom_or_hemi == "hemi":
                    raise ValueError("Unexpected state: 'hemi' variant requested on chromosome %s" % chrom)
            else:
                # sex of sample determines if it's homozygous or hemizygous
                sex = sample_id_sex[sample_id]
                if sex == "m":
                    is_male = True
                elif sex == "f":
                    is_male = False
                else:
                    raise ValueError("Sample %s has unexpected value for sex: %s" % (sample_id, sex))

                if het_or_hom_or_hemi == "hom":
                    if is_male and _is_hemizygous_segment(chrom, pos):
                        continue
                elif het_or_hom_or_hemi == "hemi":
                    if not (is_male and _is_hemizygous_segment(chrom, pos)):
                        continue

        else:
            raise ValueError("Unexpected het_or_hom_or_hemi value: " + str(het_or_hom_or_hemi))

        if DP < 10 or GQ < 20:
            continue  # skip samples that don't pass _Adj thresholds since they are not counted in the ExAC browser het/hom counts.

        if not sample_id_include_status[sample_id]:
            continue  # skip samples where include status != "YES"

        relevant_samples.append( {"sample_id": sample_id, "GQ": GQ} )

    # return samples in order from highest to lowest GQ.
    for sample in sorted(relevant_samples, key=lambda s: s["GQ"], reverse=True):
        yield sample["sample_id"]  #, sample["GQ"]

