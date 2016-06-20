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
import datetime
import getpass
import gzip
import pysam
import random
import re
import time
import traceback
from peewee import fn

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


def main(vcf_iterator, bam_output_dir, chrom=None, exit_after_minutes=None):
    """Generates HC-reassembled bams for all ExAC variants in this interval.

    Args:
        vcf_iterator: Iterator over VCF records
        bam_output_dir: Top level output dir for all bams
        chrom: (optional) genomic region to limit processing
        exit_after_minutes: (optional - integer) after this many minutes, finish processing the current variant and exit
    """

    # iterate over the VCF
    main_started_time = datetime.datetime.now()

    counters = collections.defaultdict(int)
    for chrom, pos, ref, alt_alleles, n_het_list, n_hom_list, n_hemi_list, all_genotypes_in_row in vcf_iterator:
        counters["sites"] += 1

        if exit_after_minutes:
            minutes_since_task_started = (datetime.datetime.now() - main_started_time).total_seconds()/3600
            if minutes_since_task_started > exit_after_minutes:
                logging.info("Time limit of %s minutes reached. Exiting..." % exit_after_minutes)
                break

        # iterate over alt alleles (in case this row is multi-allelic)
        for alt_allele_index, (alt, n_het, n_hom, n_hemi) in enumerate(zip(alt_alleles, n_het_list, n_hom_list, n_hemi_list)):

            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(pos, ref, alt)

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

                vr.started=1
                vr.started_time = datetime.datetime.now()
                vr.username = getpass.getuser()[0:10]
                vr.comments = str(vr.comments or "") + "_s"
                vr.save()

                if het_or_hom_or_hemi == "het":
                    n_expected_samples = n_het
                elif het_or_hom_or_hemi == "hom":
                    n_expected_samples = n_hom
                elif het_or_hom_or_hemi == "hemi":
                    n_expected_samples = n_hemi
                else:
                    raise ValueError("Unexpected value for het_or_hom_or_hemi: %s" % str(het_or_hom_or_hemi))

                if n_expected_samples == 0:
                    logging.info("%s-%s-%s-%s %s - has n_expected_samples == 0 - skipping.." % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi))
                    vr.n_expected_samples = 0
                    vr.finished = 1
                    vr.finished_time = datetime.datetime.now()
                    vr.comments = str(vr.comments or "") + "_n_expected=0"
                    vr.save()
                    continue

                # look in the genotypes vcf to get sample ids to use for readviz, sorted in order from
                # best-to-show-for-readviz to worst
                failed_sample_counter = 0
                chosen_reassembled_bams = []
                best_for_readviz_sample_ids = best_for_readviz_sample_id_iter(
                        chrom,
                        minrep_pos,
                        het_or_hom_or_hemi,
                        alt_allele_index + 1,  # + 1 because 0 means REF in the genotype objects
                        all_genotypes_in_row,
                        EXAC_SAMPLE_ID_TO_INCLUDE_STATUS,
                        EXAC_SAMPLE_ID_TO_SEX)

                # double-check that the n_hom, n_het, n_hemi
                if abs(n_expected_samples - len(best_for_readviz_sample_ids)) > 0:
                    raise ValueError("%s-%s-%s-%s %s - n_expected_samples != len(best_for_readviz_sample_ids): %s != %s" % (
                        chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi,
                        n_expected_samples, len(best_for_readviz_sample_ids)))

                for next_best_readviz_sample_id in best_for_readviz_sample_ids:
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
                        traceback.print_exc()
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
                logging.info("%s-%s-%s-%s %s - saving as finished." % (chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom_or_hemi))
                vr.n_expected_samples=min(n_expected_samples, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
                vr.n_available_samples=len(chosen_reassembled_bams)
                vr.n_failed_samples=failed_sample_counter
                vr.readviz_bam_paths="|".join(chosen_reassembled_bams)
                vr.finished=1
                vr.finished_time = datetime.datetime.now()
                vr.comments = str(vr.comments or "") + "_done_%d_of_%d" % (vr.n_available_samples, vr.n_expected_samples)
                vr.save()


    logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)
    logging.info("generate_HC_bams finished.") # at %s:%s" % (chrom, pos))


def create_iterator_from_variant_table(tabix_file, chrom=None, start_pos=None, end_pos=None):
    """Iterate over variants in the tabix file that are marked as not-yet-finished in the variant table

    Args:
        tabix_file: pysam tabix file object
        chrom: chromosome
        start_pos: integer 1-based inclusive start position of genomic region
        end_pos: integer 1-based inclusive end position of genomic region
    Returns:
        lines in the tabix file corresponding to unfinished variants
    """

    where_condition = (Variant.started == 0) & (Variant.finished == 0)
    if chrom is not None:
        where_condition = where_condition & (Variant.chrom == chrom)
    if start_pos is not None:
        where_condition = where_condition & (Variant.pos >= start_pos)
    if end_pos is not None:
        where_condition = where_condition & (Variant.pos <= end_pos)

    while True:
        # claim a variant to process. 
        with db.atomic() as txn:
            # order_by random to avoid collisions with other threads/processes
            #unprocessed_variants = Variant.select().where(where_condition).order_by(fn.Rand()).limit(1)  # order_by(fn.Rand()) queries are too slow 
            randomized_variant_num = random.randint(1, 200)
            unprocessed_variants = Variant.select().where(where_condition).limit(randomized_variant_num)  # order_by(fn.Rand()) queries are too slow 
            #logging.info("running query: " + str(unprocessed_variants.sql()))
            unprocessed_variants = list(unprocessed_variants)
            if len(unprocessed_variants) == 0:
                logging.info("Finished all variants. Exiting..")
                break

            current_variant = unprocessed_variants[-1]  # pick the last one
            logging.info("retrieving next variant..")

            query = Variant.update(started=1).where( (Variant.id == current_variant.id) & where_condition )
            logging.info("running query: " + str(query.sql()))
            rows_updated = query.execute()

        if rows_updated == 0:
            sleep_interval = random.randint(1, 15)
            logging.info("%s-%s-%s-%s - variant claimed by another task. Skipping.. (sleep for %s sec)" % (
                current_variant.chrom, current_variant.pos, current_variant.ref, current_variant.alt, sleep_interval))
            time.sleep(sleep_interval) # sleep for a random time interval to avoid constant lock contension
            continue

        # for the current_variant, get the corresponding row from the VCF
        for fields in tabix_file.fetch(str(current_variant.chrom), current_variant.pos-1, current_variant.pos+1):
            logging.info("%s-%s-%s-%s - retrieved variant " % (
                current_variant.chrom, current_variant.pos, current_variant.ref, current_variant.alt))
            if int(fields[1]) == current_variant.pos:
                yield fields
                break
        else:
            raise Exception("Variant not found in fetch results: " + str(current_variant.__dict__))


if __name__ == "__main__":
    p = configargparse.getArgumentParser()
    p.add("--bam-output-dir", help="Where to output HC-reassembled bams", default=BAM_OUTPUT_DIR)
    p.add("--chrom", help="If specified, only process this chromosome")
    p.add("--start-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int)
    p.add("--end-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int, default=10**10)

    p.add("--process-variant-table", action="store_true", help="If specified, then instead of iterating over the "
        "ExAC genotypes vcf, process unfinished variants from the SQL variant table")
    p.add("--exit-after", metavar="MINUTES", help="This many minutes after starting, finish processing "
                                                  "the current variant and then exit", default=60*3.5, type=float)

    args = p.parse_args()

    logging.info("Running with settings: ")
    for argname in filter(lambda n: not n.startswith("_"), dir(args)):
        logging.info("   %s = %s" % (argname, getattr(args, argname)))
    logging.info("\n")

    import utils.constants
    for key in dir(utils.constants):
        if not key.startswith("_"):
            logging.info("%s=%s" % (key, utils.constants.__dict__[key]))

    db = init_db()

    # initialize vcf iterator
    exac_full_vcf = EXAC_FULL_VCF_PATH
    tabix_file = pysam.TabixFile(filename=exac_full_vcf, parser=pysam.asTuple())
    last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
    exac_sample_ids = set(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS.keys())
    parse_vcf_row = create_vcf_row_parser(last_header_line, exac_sample_ids)

    if args.process_variant_table:
        logging.info("Processing remaining variants in variant table")

        vcf_row_iterator = create_iterator_from_variant_table(tabix_file, args.chrom, args.start_pos, args.end_pos)
    elif args.chrom:
        logging.info("Processing variants in %s:%s-%s in %s" % (args.chrom, args.start_pos, args.end_pos, exac_full_vcf))

        if args.start_pos:
            args.start_pos = args.start_pos - 1  # because start_pos is 1-based inclusive and fetch(..) doesn't include the start_pos

        vcf_row_iterator = tabix_file.fetch(args.chrom, args.start_pos, args.end_pos)
    else:
        logging.info("Processing all variants in %s" % exac_full_vcf)
        vcf_row_iterator = pysam.tabix_iterator(gzip.open(exac_full_vcf), parser=pysam.asTuple())

    #if args.start_pos:
    #    vcf_row_iterator = (row for row in vcf_row_iterator if int(row[1]) >= args.start_pos)
    #if args.end_pos:
    #    vcf_row_iterator = (row for row in vcf_row_iterator if int(row[1]) <= args.end_pos)

    vcf_iterator = (parse_vcf_row(row) for row in vcf_row_iterator)

    main(vcf_iterator=vcf_iterator,
         bam_output_dir=args.bam_output_dir,
         chrom=args.chrom,
         exit_after_minutes=args.exit_after
    )
