"""
Utility methods for choosing bam samples to display for a particular ExAC variant.
"""

# how many samples to show per het or hom-alt variant in the exac browser.
MAX_SAMPLES_TO_SHOW_PER_VARIANT = 5

# in ExAC v3 the max ref allele is 325bp long, and max alt allele is 365bp long, so use 500 to be safe.
MAX_ALLELE_SIZE = 500  # nucleotides


import os
import re
import time


def get_bam_path(sample_id, sample_id_to_bam_path):
    """This function tries to find the bam path for the given sample id.

    Args:
      sample_id: vcf sample id
      sample_id_to_bam_path: based on exac info table
    Return: The bam path if the bam file exists.
    Raise: IOError if the bam doesn't exist.
    """

    # work-arounds for relocated .bams based on @birndle's igv_spot_checking script
    bam_path = sample_id_to_bam_path[sample_id]
    if "/cga/pancan2/picard_bams/ext_tcga" in bam_path:
        bam_path = bam_path.replace("/cga/pancan2/", "/cga/fh/cga_pancan2/")

    if "CONT_" in bam_path:
        bam_path = bam_path.replace("CONT_", "CONT")

    bam_path = re.sub("/v[0-9]{1,2}/", "/current/", bam_path)  # get latest version of the bam

    return bam_path



def choose_samples(het_or_hom, alt_allele_index, genotypes, sample_id_include_status, sample_id_to_bam_path, sample_id_to_gvcf_path):
    """Contains heuristics for choosing which samples to display for a given
    variant (chrom, pos, ref, alt) in it's het or hom-alt state.

    Args:
        het_or_hom: Either "het" or "hom" to indicate whether to choose het or hom-alt samples.
        alt_allele_index: 1-based index of the alt allele that we're choosing samples for (bewteen 1 and n = total number of alt alleles at this site)
        genotypes: a dictionary that maps each sample_id to a 4-tuple: (gt_ref, gt_alt, GQ, DP)
            where gt_ref and gt_alt are integers between 0, and num alt alleles.
        sample_id_include_status: dictionary mapping sample_id to either
            True or False depending on the last column of the exac info table
        sample_id_to_bam_path: bam path from the ExAC info table

    Return:
        Three lists containing file paths: chosen_sample_bams, chosen_sample_gvcfs, missing_bams
        The list length of chosen_sample_bams and chosen_sample_gvcfs will be <=
        MAX_SAMPLES_TO_SHOW_PER_VARIANT. missing_bams will be non-empty if there
        any samples that would have been chosen, except their bam is missing on the filesystem
    """
    assert het_or_hom in ["het", "hom"], "Unexpected het_or_hom arg: %s" % het_or_hom

    # filter down from all samples to just the samples that have the desired genotype and have include=YES
    relevant_samples = []  # a list of dicts
    for sample_id, (gt_ref, gt_alt, GQ, DP) in genotypes.items():
        if gt_ref is None and gt_alt is None:
            continue

        if het_or_hom == "het":
            if gt_ref == gt_alt:
                continue  # skip unless heterozygous
            if gt_ref != alt_allele_index and gt_alt != alt_allele_index:
                continue # check both gt_ref and gt_alt here to handle het genotypes that are alt-alt  (eg. 1/2)
        elif het_or_hom == "hom":
            if gt_ref != gt_alt:
                continue  # skip unless homozygous
            if gt_alt != alt_allele_index:
                continue  # skip unless homozygous for the specific alt allele we're looking for (this matters for multiallelics)
        else:
            raise ValueError("Unexpected het_or_hom value: " + het_or_hom)

        if DP < 10 or GQ < 20:
            continue  # ignore samples that don't pass _Adj thresholds since they are not counted in the ExAC browser het/hom counts.
        if sample_id_include_status[sample_id]:
            relevant_samples.append( {"sample_id": sample_id, "GQ": GQ} )


    # get up to MAX_SAMPLES_TO_SHOW_PER_VARIANT samples with the highest GQ.
    # skip samples whose bams aren't found on disk.
    relevant_samples.sort(key=lambda s: s["GQ"], reverse=True)

    # figure out list of bam paths for samples to display
    chosen_sample_bams = []
    chosen_sample_gvcfs = []
    missing_bams = []
    while relevant_samples and len(chosen_sample_bams) < MAX_SAMPLES_TO_SHOW_PER_VARIANT:
        # retrieve the sample with the next-highest GQ, and remove it from the
        # list so its not considered again
        max_GQ_sample = relevant_samples[0]
        del relevant_samples[0]

        # convert sample id to bam path
        sample_id = max_GQ_sample["sample_id"]
        bam_path = get_bam_path(sample_id, sample_id_to_bam_path)

        # check if this bam exists
        for retry_counter in range(3):
            if os.access(bam_path, os.R_OK):
                chosen_sample_bams.append(bam_path)

                gvcf_path = sample_id_to_gvcf_path[sample_id]
                chosen_sample_gvcfs.append(gvcf_path)
                break
            time.sleep(1)
        else:
            missing_bams.append(bam_path)

    return chosen_sample_bams, chosen_sample_gvcfs, missing_bams
