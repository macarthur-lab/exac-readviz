"""
Utility methods for choosing bam samples to display for a particular ExAC variant.
"""
import logging
from collections import OrderedDict

def _is_hemizygous_segment(chrom, pos):
    """Utility method that takes a chromosome (eg. '1', '2', 'X'..) and a position integer and returns True if the
     it's on the X or Y chromosome and not in the PAR region.
    """
    return (chrom == 'X' and not ((60001 <= pos <= 2699520) or (154931044 <= pos <= 155260560))) or (chrom == 'Y')


def _is_male(sex):
    """Utility method that returns True if sex == "m", returns False if sex == "f", and raises an error otherwise. """
    if sex == "m":
        return True
    elif sex == "f":
        return False
    else:
        raise ValueError("Unexpected value for sex arg: %s" % sex)


def _is_hemizygous(chrom, pos, sex):
    """Utility method that takes chrom (eg. '1', 'X', 'Y', etc.), position as an integer, and sex as either "m" or "f"
    and returns True if this position should be considered hemizygous"""
    return chrom in ('X', 'Y') and _is_male(sex) and _is_hemizygous_segment(chrom, pos)


def best_for_readviz_sample_id_iter(chrom, pos, het_or_hom_or_hemi, alt_allele_index, genotypes, sample_id_include_status, sample_id_sex):
    """Implements heuristics for choosing which samples are best to display for
    a given variant (chrom, pos, ref, alt) in it's het or hom-alt state.

    Args:
        chrom: chromosome string "chr1", "chrX", etc.
        pos:  variant position as an integer (minified)
        het_or_hom_or_hemi: Either "het" or "hom" or "hemi" to indicate whether to choose het, hom-alt, or hemi samples.
        alt_allele_index: 1-based index of the alt allele that we're choosing
            samples for (bewteen 1 and n = total number of alt alleles at this site)
        genotypes: a dictionary that maps each VCF sample_id to a 4-tuple: (gt_ref, gt_alt, AD, DP, GQ)
            where gt_ref and gt_alt are integers between 0 and num alt alleles.
        sample_id_include_status: dictionary mapping each VCF sample_id to either
            True or False depending on the last column of the exac info table
        sample_id_sex: dictionary mapping each VCF sample_id to "m" or "f"

    Return:
        Iterator over sample ids in order - from those that should be first to
        be displayed to those that should be last.
    """
    assert het_or_hom_or_hemi in ["het", "hom", "hemi"], "Unexpected het_or_hom_or_hemi arg: %s" % het_or_hom_or_hemi

    counter = OrderedDict()
    
    # filter all samples down to just samples that have the desired genotype and have include=YES
    relevant_samples = []  # a list of dicts
    for sample_id, (gt_ref, gt_alt, AD, DP, GQ) in genotypes.items():
        counter['total samples'] = counter.get('total samples', 0) + 1
        if gt_ref is None and gt_alt is None:
            continue
        counter['sample genotype is not ./.'] = counter.get('sample genotype is not ./.', 0) + 1

        if het_or_hom_or_hemi == "het":
            if gt_ref == gt_alt:
                continue  # skip if homozygous
            if gt_ref != alt_allele_index and gt_alt != alt_allele_index:
                continue # skip if neither allele matches the specific alt allele we're looking for (eg. 1/3)
            if _is_hemizygous(chrom, pos, sample_id_sex[sample_id]):
                continue
        elif het_or_hom_or_hemi == "hom":
            if gt_ref != gt_alt:
                continue  # skip unless homozygous
            if gt_alt != alt_allele_index:
                continue  # skip unless homozygous for the specific alt allele we're looking for (this matters for multiallelics)
            if _is_hemizygous(chrom, pos, sample_id_sex[sample_id]):
                continue
        elif het_or_hom_or_hemi == "hemi" and chrom in ('X', 'Y'):
            if gt_ref != alt_allele_index and gt_alt != alt_allele_index:
                continue # skip if neither allele matches the specific alt allele we're looking for (eg. 1/3)
            counter['correct allele'] = counter.get('correct allele', 0) + 1

            if not _is_hemizygous(chrom, pos, sample_id_sex[sample_id]):
                continue
            counter['sample sex == male and chrom, pos in hemizigous_segment'] = counter.get('sample sex == male and chrom, pos in hemizigous_segment', 0) + 1

            # check whether the allele with the most reads is the ref allele, in which case it's homozygous reference
            if AD is None:
                continue
            if gt_ref != gt_alt:
                if len(AD) > 2:
                    # handle multiallelics using Monkol's AC_Hemi code
                    max_ad = AD[gt_ref]
                    hemi_allele_i = gt_ref
                    for i in range(0, len(AD)):
                        if AD[i] > max_ad:
                            hemi_allele_i = i
                    if hemi_allele_i != gt_ref:
                        hemi_allele_i = gt_alt
                    if hemi_allele_i != gt_alt:
                        continue
                elif AD[0] >= AD[1]:
                    continue

            counter['AD alt > AD ref'] = counter.get('AD alt > AD ref', 0) + 1
            logging.info("%s  %s/%s:%s:%s:%s " % (sample_id, gt_ref, gt_alt, ",".join(map(str, AD)), DP, GQ))
        else:
            raise ValueError("Unexpected het_or_hom_or_hemi value: " + str(het_or_hom_or_hemi))

        if DP < 10 or GQ < 20:
            continue  # skip samples that don't pass _Adj thresholds since they are not counted in the ExAC browser het/hom counts.
        counter["sample passes DP>=10 and GQ>=20"] = counter.get('sample passes DP>=10 and GQ>=20', 0) + 1

        if not sample_id_include_status[sample_id]:
            continue  # skip samples where include status != "YES"
        counter["INCLUDE status == YES"] = counter.get('INCLUDE status == YES', 0) + 1

        relevant_samples.append( {"sample_id": sample_id, "GQ": GQ} )

    if het_or_hom_or_hemi == "hemi":
        for k in counter:
            logging.info("   --- hemi counts: %s: %s" % (counter[k], k))

    # return sample ids in order from highest to lowest GQ.
    return [sample["sample_id"] for sample in sorted(relevant_samples, key=lambda s: s["GQ"], reverse=True)]


