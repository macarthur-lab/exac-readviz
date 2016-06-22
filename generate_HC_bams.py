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
import signal
import time
import traceback
from peewee import fn

from utils.database import init_db, Variant
from utils.choose_samples import best_for_readviz_sample_id_iter
from utils.constants import MAX_SAMPLES_TO_SHOW_PER_VARIANT, EXAC_FULL_VCF_PATH, BAM_OUTPUT_DIR
from utils.exac_info_table import EXAC_SAMPLE_ID_TO_GVCF_PATH, EXAC_SAMPLE_ID_TO_INCLUDE_STATUS, EXAC_SAMPLE_ID_TO_SEX
from utils.exac_info_table import lookup_original_bam_path
from utils.exac_vcf import create_vcf_row_parser, create_variant_iterator_from_vcf
from utils.haplotype_caller import run_haplotype_caller


import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p')


CTRL_C_SIGNAL = False
def signal_handler(signal, frame):
    global CTRL_C_SIGNAL
    CTRL_C_SIGNAL = True
    logging.info("Ctrl-C pressed")
    
signal.signal(signal.SIGINT, signal_handler)


def main(variant_iterator, bam_output_dir, exit_after_minutes=None):
    """Generates HC-reassembled bams for all ExAC variants in this interval.

    Args:
        variant_iterator: Iterator that returns exac_vcf.VariantGT objects
        bam_output_dir: Top level output dir for all bams
        exit_after_minutes: (optional - integer) after this many minutes, finish processing the current variant and exit
    """
    
    # iterate over the VCF
    main_started_time = datetime.datetime.now()

    counters = collections.defaultdict(int)
    for chrom, pos, ref, alt, alt_allele_index, het_or_hom_or_hemi, n_expected_samples, all_genotypes_in_row in variant_iterator:
        # print some stats
        logging.info("-----")

        counters["all_variants_gt"] += 1
        if counters["all_variants_gt"] % 100 == 0:
            logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)
            logging.info("-----")

        if exit_after_minutes:
            minutes_since_task_started = (datetime.datetime.now() - main_started_time).total_seconds()/3600
            if minutes_since_task_started > exit_after_minutes:
                logging.info("Time limit of %s minutes reached. Exiting..." % exit_after_minutes)
                break
        if CTRL_C_SIGNAL:
            logging.info("Interrupted. Exiting...")
            break

        # check if allele has been processed already (skip if yes)
        vr, created = Variant.get_or_create(chrom=chrom, pos=pos, ref=ref, alt=alt, het_or_hom_or_hemi=het_or_hom_or_hemi)

        if vr.finished:
            logging.info("%s-%s-%s-%s %s - already done (%s out of %s available) - skipping.." % (
                chrom, pos, ref, alt, het_or_hom_or_hemi, vr.n_available_samples,
                vr.n_expected_samples))
            counters[het_or_hom_or_hemi+"_variants_gt_already_done"] += 1
            continue

        vr.variant_id = "%s-%s-%s-%s" % (chrom, pos, ref, alt)
        vr.started = 1
        vr.started_time = datetime.datetime.now()
        vr.username = getpass.getuser()[0:10]
        vr.comments = str(vr.comments or "") + "_s"
        vr.save()

        if n_expected_samples == 0:
            logging.info("%s-%s-%s-%s %s - has n_expected_samples == 0 - skipping.." % (chrom, pos, ref, alt, het_or_hom_or_hemi))
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
                pos,
                het_or_hom_or_hemi,
                alt_allele_index + 1,  # + 1 because 0 means REF in the genotype objects
                all_genotypes_in_row,
                EXAC_SAMPLE_ID_TO_INCLUDE_STATUS,
                EXAC_SAMPLE_ID_TO_SEX)

        # double-check that the n_hom, n_het, n_hemi
        if abs(n_expected_samples - len(best_for_readviz_sample_ids)) > 0:
            raise ValueError("%s-%s-%s-%s %s - n_expected_samples != len(best_for_readviz_sample_ids): %s != %s" % (
                chrom, pos, ref, alt, het_or_hom_or_hemi, n_expected_samples, len(best_for_readviz_sample_ids)))

        for next_best_readviz_sample_id in best_for_readviz_sample_ids:
            try:
                original_bam_path = lookup_original_bam_path(next_best_readviz_sample_id)
                original_gvcf_path = EXAC_SAMPLE_ID_TO_GVCF_PATH[next_best_readviz_sample_id]

                succeeded, reassembled_bam_path = run_haplotype_caller(
                    chrom,
                    pos,
                    ref,
                    alt,
                    het_or_hom_or_hemi,
                    original_bam_path,
                    original_gvcf_path,
                    bam_output_dir,
                    next_best_readviz_sample_id,
                    sample_i=len(chosen_reassembled_bams))
            except Exception as e:
                logging.error("%s-%s-%s-%s %s - error in run_haplotype_caller: %s" % (chrom, pos, ref, alt, het_or_hom_or_hemi, e))
                traceback.print_exc()
                succeeded = False

            if succeeded:
                chosen_reassembled_bams.append(reassembled_bam_path)
                if len(chosen_reassembled_bams) >= n_expected_samples:
                    logging.info("%s-%s-%s-%s %s - all n_expected_samples == %d now found." % (chrom, pos, ref, alt, het_or_hom_or_hemi, n_expected_samples))
                    break
                if len(chosen_reassembled_bams) >= MAX_SAMPLES_TO_SHOW_PER_VARIANT:
                    logging.info("%s-%s-%s-%s %s - all %d samples now found." % (chrom, pos, ref, alt, het_or_hom_or_hemi, MAX_SAMPLES_TO_SHOW_PER_VARIANT))
                    break
            else:
                failed_sample_counter += 1

        # save variant record
        logging.info("%s-%s-%s-%s %s - saving as finished." % (chrom, pos, ref, alt, het_or_hom_or_hemi))
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


def create_variant_record_iterator(chrom=None, start_pos=None, end_pos=None):
    """Iterate over variant records that are marked as not-yet-finished in the variant table

    Args:
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
        # get the next variant to process
        randomized_variant_num = random.randint(1, 200)
        unprocessed_variants = Variant.select().where(where_condition).limit(randomized_variant_num)  # order_by(fn.Rand()).limit(1) queries are too slow
        unprocessed_variants_list = list(unprocessed_variants)
        if len(unprocessed_variants_list) == 0:
            logging.info("Finished all variants. Exiting..")
            break

        #sql, params = unprocessed_variants.sql()
        #logging.info("query: %s \n rows retreived %s" % ( (sql % tuple(params)), len(unprocessed_variants_list)))

        current_variant = unprocessed_variants_list[-1]  # pick the last one
        logging.info("retrieving next variant id = %s: %s-%s-%s-%s %s" % (current_variant.id,
            current_variant.chrom, current_variant.pos, current_variant.ref, current_variant.alt, current_variant.het_or_hom_or_hemi))

        with db.atomic():
            query = Variant.update(started=1).where( (Variant.id == current_variant.id) & where_condition )
            rows_updated = query.execute()

        #sql, params = query.sql()
        #logging.info("query: %s \n rows updated: %s" % ( (sql % tuple(params)), rows_updated))

        if rows_updated == 0:
            sleep_interval = random.randint(1, 15)
            logging.info("%s-%s-%s-%s - variant claimed by another task. Skipping.. (sleep for %s sec)" % (
                current_variant.chrom, current_variant.pos, current_variant.ref, current_variant.alt, sleep_interval))
            time.sleep(sleep_interval) # sleep for a random time interval to avoid constant lock contension
            continue

        yield current_variant


def create_variant_iterator_from_variant_records(variant_record_iterator, tabix_file, parse_vcf_row):
    """Takes a variant_record_iterator and yields VariantGT objects based on corresponding rows parsed from the VCF"""

    for variant_record in variant_record_iterator:
        # for the current_variant, read the corresponding row from the VCF
        current_variant_record = current_variant_vcf_fields = None
        for fields in tabix_file.fetch(str(variant_record.chrom), variant_record.pos-1, variant_record.pos+1):
            logging.info("%s-%s-%s-%s %s - retrieved variant " % (
                variant_record.chrom, variant_record.pos, variant_record.ref, variant_record.alt, variant_record.het_or_hom_or_hemi))
            if int(fields[1]) == variant_record.pos:
                current_variant_vcf_fields = fields
                current_variant_record = variant_record
                break
        else:
            raise Exception("Variant not found in fetch results: %s" % str(variant_record.__dict__))

        # parse this row and, since it may represent multiple variants, yield one or more VariantGT objects from it
        for variant in parse_vcf_row(current_variant_vcf_fields, het_or_hom_or_hemi=current_variant_record.het_or_hom_or_hemi):
            yield variant




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

    if args.start_pos:
        args.start_pos = args.start_pos - 1  # because start_pos is 1-based inclusive and fetch(..) doesn't include the start_pos

    db = init_db()

    # initialize vcf row parser
    exac_full_vcf = EXAC_FULL_VCF_PATH
    tabix_file = pysam.TabixFile(filename=exac_full_vcf, parser=pysam.asTuple())
    last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
    exac_sample_ids = set(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS.keys())
    parse_vcf_row = create_vcf_row_parser(last_header_line, exac_sample_ids)

    # create variant iterator
    if args.process_variant_table:
        # from variant table
        logging.info("Processing remaining variants in variant table")
        variant_record_iterator = create_variant_record_iterator(chrom=args.chrom, start_pos=args.start_pos, end_pos=args.end_pos)

        variant_iterator = create_variant_iterator_from_variant_records(variant_record_iterator, tabix_file, parse_vcf_row)
    else:
        # from vcf
        if args.chrom:
            logging.info("Processing variants in %s:%s-%s in %s" % (args.chrom, args.start_pos, args.end_pos, exac_full_vcf))
            vcf_row_iterator = tabix_file.fetch(args.chrom, args.start_pos, args.end_pos)
        else:
            logging.info("Processing all variants in %s" % exac_full_vcf)
            vcf_row_iterator = pysam.tabix_iterator(gzip.open(exac_full_vcf), parser=pysam.asTuple())

        variant_iterator = create_variant_iterator_from_vcf(vcf_row_iterator, parse_vcf_row)

    # process the variants
    main(variant_iterator=variant_iterator,
         bam_output_dir=args.bam_output_dir,
         exit_after_minutes=args.exit_after
    )
