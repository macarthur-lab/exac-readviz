"""
This script generates HC-reassembled bams. For each variant, it
selects which samples will be displayed and then runs HaplotypeCaller on those
samples. An unlimited number of instances of this script can be
run in parallel (for example in an array job) to speed up execution.
"""


import configargparse
configargparse.initArgumentParser(
    default_config_files=["~/.generate_HC_bams_config"],
    formatter_class=configargparse.ArgumentDefaultsHelpFormatter)

import collections
import gzip
import pysam
import re

from utils.database import init_db, Variant
from utils.choose_samples import best_for_readviz_sample_id_iter
from utils.constants import MAX_SAMPLES_TO_SHOW_PER_VARIANT, EXAC_FULL_VCF_PATH, BAM_OUTPUT_DIR
from utils.exac_info_table import EXAC_SAMPLE_ID_TO_BAM_PATH, EXAC_SAMPLE_ID_TO_GVCF_PATH, \
    EXAC_SAMPLE_ID_TO_INCLUDE_STATUS, EXAC_SAMPLE_ID_TO_SEX
from utils.exac_vcf import create_vcf_row_parser
from utils.minimal_representation import get_minimal_representation
from utils.haplotype_caller import run_haplotype_caller


import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p')


def lookup_original_bam_path(sample_id):
    """Look up the bam path for the given sample id.

    Args:
      sample_id: vcf sample id
    Return:
      The original bam path
    """

    # work-arounds for relocated .bams based on @birndle's igv_spot_checking script
    bam_path = EXAC_SAMPLE_ID_TO_BAM_PATH[sample_id]
    if "/cga/pancan2/picard_bams/ext_tcga" in bam_path:
        bam_path = bam_path.replace("/cga/pancan2/", "/cga/fh/cga_pancan2/")

    if "CONT_" in bam_path:
        bam_path = bam_path.replace("CONT_", "CONT")

    bam_path = re.sub("/v[0-9]{1,2}/", "/current/", bam_path)  # get latest version of the bam

    return bam_path


def main(exac_full_vcf, bam_output_dir, chrom=None, start_pos=None, end_pos=10**10):
    """Generates HC-reassembled bams for all ExAC variants in this inteveral.

    Args:
        exac_full_vcf: The path of the full vcf
        bam_output_dir: Top level output dir for all bams
        chrom: optional genomic region to limit processing
        start_pos: integer 1-based inclusive start position of genomic region
        end_pos: integer 1-based inclusive end position of genomic region
    """

    # parse the VCF
    tabix_file = pysam.TabixFile(filename=exac_full_vcf, parser=pysam.asTuple())
    last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
    if chrom:
        logging.info("Processing variants in %s:%s-%s in %s" % (
            chrom, start_pos, end_pos, exac_full_vcf))

        if start_pos:
            start_pos = start_pos - 1  # because start_pos is 1-based inclusive and fetch(..) doesn't include the start_pos

        vcf_iterator = tabix_file.fetch(chrom, start_pos, end_pos)
    else:
        logging.info("Processing all variants in %s" % exac_full_vcf)
        vcf_iterator = pysam.tabix_iterator(gzip.open(exac_full_vcf), parser=pysam.asTuple())

    exac_sample_ids = set(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS.keys())
    parse_vcf_row = create_vcf_row_parser(last_header_line, exac_sample_ids)

    counters = collections.defaultdict(int)
    for row in vcf_iterator:
        counters["sites"] += 1

        chrom, pos, ref, alt_alleles, n_het_list, n_hom_list, n_hemi_list, all_genotypes_in_row = parse_vcf_row(row)

        # iterate over alt alleles (in case this row is multi-allelic)
        for alt_allele_index, (alt, n_het, n_hom, n_hemi) in enumerate(zip(alt_alleles, n_het_list, n_hom_list, n_hemi_list)):

            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(
                pos, ref, alt)

            # print some stats
            logging.info("-----")

            counters["all_alleles"] += 1
            if counters["all_alleles"] % 100 == 0:
                logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)

            # choose het, hom-alt, and hemizygous samples to display for this allele
            possible_genotypes = ["het", "hom", "hemi"] if chrom in ('X', 'Y') else ["het", "hom"]
            for het_or_hom_or_hemi in possible_genotypes:

                # check if allele has been processed already (skip if yes)
                vr, created = Variant.get_or_create(
                    chrom=chrom,
                    pos=minrep_pos,
                    ref=minrep_ref,
                    alt=minrep_alt,
                    het_or_hom_or_hemi = het_or_hom_or_hemi)

                if vr.finished:
                    logging.info("%s-%s-%s-%s %s - already done (%s out of %s available) - skipping.." % (
                        chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi, vr.n_available_samples,
                        vr.n_expected_samples))
                    counters[het_or_hom_or_hemi+"_alleles_already_done"] += 1
                    continue

                n_expected_samples = n_het if het_or_hom_or_hemi == "het" else n_hom
                if n_expected_samples == 0:
                    logging.info("%s-%s-%s-%s %s - has n_expected_samples == 0 - skipping.." % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi))
                    vr.n_expected_samples=0
                    vr.finished=1
                    vr.save()
                    continue

                # iterate over samples in order from best-to-show-for-readviz to worst
                failed_sample_counter = 0
                chosen_reassembled_bams = []
                for next_best_readviz_sample_id in best_for_readviz_sample_id_iter(
                    chrom,
                    minrep_pos,
                    het_or_hom_or_hemi,
                    alt_allele_index + 1,  # + 1 because 0 means REF in the genotype objects
                    all_genotypes_in_row,
                    EXAC_SAMPLE_ID_TO_INCLUDE_STATUS,
                    EXAC_SAMPLE_ID_TO_SEX):

                    try:
                        original_bam_path = lookup_original_bam_path(next_best_readviz_sample_id)
                        original_gvcf_path = EXAC_SAMPLE_ID_TO_GVCF_PATH[next_best_readviz_sample_id]

                        succeeded, reassembled_bam_path = run_haplotype_caller(
                            chrom,
                            minrep_pos,
                            minrep_ref,
                            minrep_alt,
                            het_or_hom_or_hemi,
                            original_bam_path,
                            original_gvcf_path,
                            bam_output_dir,
                            next_best_readviz_sample_id,
                            sample_i=len(chosen_reassembled_bams))
                    except Exception as e:
                        logging.error("%s-%s-%s-%s %s - error in run_haplotype_caller: %s" % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi, e))
                        succeeded = False

                    if succeeded:
                        chosen_reassembled_bams.append(reassembled_bam_path)
                        if len(chosen_reassembled_bams) >= n_expected_samples:
                            logging.info("%s-%s-%s-%s %s - all n_expected_samples == %d now found." % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi, n_expected_samples))
                            break
                        if len(chosen_reassembled_bams) >= MAX_SAMPLES_TO_SHOW_PER_VARIANT:
                            logging.info("%s-%s-%s-%s %s - all %d samples now found." % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi, MAX_SAMPLES_TO_SHOW_PER_VARIANT))
                            break
                    else:
                        failed_sample_counter += 1

                # save variant record
                vr.n_expected_samples=min(n_expected_samples, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
                vr.n_available_samples=len(chosen_reassembled_bams)
                vr.n_failed_samples=failed_sample_counter
                vr.readviz_bam_paths="|".join(chosen_reassembled_bams)
                vr.finished=1
                vr.save()


    logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)
    logging.info("generate_HC_bams finished: %s:%s-%s" % (chrom, start_pos, end_pos))

if __name__ == "__main__":
    p = configargparse.getArgumentParser()
    p.add("--bam-output-dir", help="Where to output HC-reassembled bams", default=BAM_OUTPUT_DIR)
    p.add("--chrom", help="If specified, only process this chromosome")
    p.add("--start-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int)
    p.add("--end-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int)
    args = p.parse_args()

    logging.info("Running with settings: ")
    for argname in filter(lambda n: not n.startswith("_"), dir(args)):
        logging.info("   %s = %s" % (argname, getattr(args, argname)))
    logging.info("\n")

    import utils.constants
    for key in dir(utils.constants):
        if not key.startswith("_"):
            logging.info("%s=%s" % (key, utils.constants.__dict__[key]))

    init_db()

    main(exac_full_vcf=EXAC_FULL_VCF_PATH,
         bam_output_dir=args.bam_output_dir,
         chrom=args.chrom,
         start_pos=args.start_pos,
         end_pos=args.end_pos,
    )
