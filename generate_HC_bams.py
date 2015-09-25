"""
This script generates HC-reassembled bams. For each variant, it
selects which samples will be displayed and then runs HaplotypeCaller on those
samples. An unlimited number of instances of this script can be
run in parallel (for example in an array job) to speed up execution.
"""


import configargparse
configargparse.initArgumentParser(default_config_files=["~/.generate_HC_bams_config"],
                                  formatter_class=configargparse.ArgumentDefaultsHelpFormatter)

import collections
import datetime
import gzip
import logging
import os
import pysam
import re
import subprocess

from utils.database import init_db, Variant
from utils.choose_samples import best_for_readviz_sample_id_iter
from utils.constants import MAX_SAMPLES_TO_SHOW_PER_VARIANT, EXAC_FULL_VCF_PATH
from utils.exac_info_table import EXAC_SAMPLE_ID_TO_BAM_PATH, EXAC_SAMPLE_ID_TO_GVCF_PATH, EXAC_SAMPLE_ID_TO_INCLUDE_STATUS
from utils.exac_vcf import create_vcf_row_parser
from utils.minimal_representation import get_minimal_representation
from utils.haplotype_caller import run_haplotype_caller


import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


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

    # parse the VCF
    tabix_file = pysam.TabixFile(filename=exac_full_vcf, parser=pysam.asTuple())
    last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
    if chrom:
        logging.info("Processing variants in %s:%s-%s in %s" % (
            chrom, start_pos, end_pos, exac_full_vcf))
        vcf_iterator = tabix_file.fetch(chrom, start_pos, end_pos)
    else:
        logging.info("Processing all variants in %s" % exac_full_vcf)
        vcf_iterator = pysam.tabix_iterator(gzip.open(exac_full_vcf), parser=pysam.asTuple())

    exac_sample_ids = set(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS.keys())
    parse_vcf_row = create_vcf_row_parser(last_header_line, exac_sample_ids)

    counters = collections.defaultdict(int)
    for row in vcf_iterator:
        counters["sites"] += 1

        chrom, pos, ref, alt_alleles, n_het_list, n_hom_list, all_genotypes_in_row = parse_vcf_row(row)

        # iterate over alt alleles (in case this row is multi-allelic)
        for alt_allele_index, (alt, n_het, n_hom) in enumerate(zip(alt_alleles, n_het_list, n_hom_list)):

            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(
                pos, ref, alt)

            # print some stats
            counters["all_alleles"] += 1
            if counters["all_alleles"] % 100 == 0:
                logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)

            # choose het and hom-alt samples to display for this allele
            for het_or_hom in ["het", "hom"]:
                # check if allele has been processed already (skip if yes)
                vr, created = Variant.get_or_create(
                    chrom=chrom,
                    minrep_pos=minrep_pos,
                    minrep_ref=minrep_ref,
                    minrep_alt=minrep_alt,
                    het_or_hom=het_or_hom)

                if vr.finished:
                    logging.info("%s-%s-%s-%s %s - already done" % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom))
                    counters[het_or_hom+"_alleles_already_done"] += 1
                    continue

                # iterate over samples in order from best-to-show-for-readviz to worst
                chosen_reassembled_bams = []
                for next_best_readviz_sample_id in best_for_readviz_sample_id_iter(
                            het_or_hom,
                            alt_allele_index + 1,  # + 1 because 0 means REF in the genotype objects
                            all_genotypes_in_row,
                            EXAC_SAMPLE_ID_TO_INCLUDE_STATUS):

                    original_bam_path = lookup_original_bam_path(next_best_readviz_sample_id)
                    original_gvcf_path = EXAC_SAMPLE_ID_TO_GVCF_PATH[next_best_readviz_sample_id]

                    succeeded, reassembled_bam_path = run_haplotype_caller(
                        chrom,
                        minrep_pos,
                        minrep_ref,
                        minrep_alt,
                        het_or_hom,
                        original_bam_path,
                        original_gvcf_path,
                        bam_output_dir,
                        next_best_readviz_sample_id,
                        sample_i=len(chosen_reassembled_bams))

                    if succeeded:
                        chosen_reassembled_bams.append(reassembled_bam_path)
                        if len(chosen_reassembled_bams) >= MAX_SAMPLES_TO_SHOW_PER_VARIANT:
                            break

                # save variant record
                vr.n_expected_samples=n_het if het_or_hom == "het" else n_hom
                vr.n_available_samples=len(chosen_reassembled_bams)
                vr.readviz_bam_paths="|".join(chosen_reassembled_bams)
                vr.finished=1
                vr.save()


    logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)
    logging.info("Finished processing interval: %s:%s-%s" % (chrom, start_pos, end_pos))

if __name__ == "__main__":
    p = configargparse.getArgumentParser()
    p.add("--bam-output-dir", help="Where to output HC-reassembled bams", default='/broad/hptmp/exac_readviz_backend/')
    p.add("--chrom", help="If specified, only process this chromosome")
    p.add("--start-pos", help="If specified, only process region in this interval", type=int)
    p.add("--end-pos", help="If specified, only process region in this interval", type=int)
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
