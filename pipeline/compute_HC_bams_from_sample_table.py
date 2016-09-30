"""
This script generates HC-reassembled bams by processing records in the Sample table
"""

import logging
logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')

logging.info("compute_HC_bams_from_sample_table - running")

import configargparse

from utils.exac_info_table import lookup_original_bam_path

configargparse.initArgumentParser(
        default_config_files=["~/.generate_HC_bams_config"],
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)

import collections
import datetime
import pysam
import random
import signal
import time
import traceback

logging.info("compute_HC_bams_from_sample_table - done with imports - #1")

from utils.database import init_db, Sample, _readviz_db
logging.info("compute_HC_bams_from_sample_table - done with imports - #2")

from utils.constants import BAM_OUTPUT_DIR, MAX_SAMPLES_TO_SHOW_PER_VARIANT, BACKUP_SAMPLES_IN_CASE_OF_ERRORS, HUMAN_GENOME_FASTA_PATH
PYSAM_FASTA = pysam.FastaFile(HUMAN_GENOME_FASTA_PATH)

logging.info("compute_HC_bams_from_sample_table - done with imports - #3")

from utils.haplotype_caller import run_haplotype_caller
logging.info("compute_HC_bams_from_sample_table - done with imports - #4")

from utils.normalize import normalize

CTRL_C_SIGNAL = False
def signal_handler(signal, frame):
    global CTRL_C_SIGNAL
    CTRL_C_SIGNAL = True
    logging.info("Ctrl-C pressed")

signal.signal(signal.SIGINT, signal_handler)


def main(sample_iterator, bam_output_dir, exit_after_minutes=None):
    """Generates HC-reassembled bams.

    Args:
        sample_iterator: Iterator that returns Sample records.
        bam_output_dir: Top level output dir for all bams
        exit_after_minutes: (optional - integer) after this many minutes, finish processing the current sample and exit
    """

    # iterate over the Samples
    main_started_time = datetime.datetime.now()

    counters = collections.defaultdict(int)
    for sr in sample_iterator:
        # print some stats
        logging.info("-----")

        counters["all_samples_gt"] += 1
        if counters["all_samples_gt"] % 100 == 0:
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

        # skip if sample has been processed already -- this should never happen
        if sr.finished:
            logging.info("%s-%s-%s-%s %s - already done - skipping.." % (
                sr.chrom, sr.pos, sr.ref, sr.alt, sr.het_or_hom_or_hemi))
            counters[sr.het_or_hom_or_hemi+"_sample_already_done"] += 1
            continue

        # compute sample_i if it hasn't been already - this doesn't work that well because race condition may (on rare occasion) cause 2 records to have the same sample_i
        if sr.sample_i is None:
            available_sample_ids = [
                s.id for s in Sample.select().where(
                    (Sample.chrom == sr.chrom) &
                    (Sample.pos == sr.pos) &
                    (Sample.ref == sr.ref) &
                    (Sample.alt == sr.alt) &
                    (Sample.het_or_hom_or_hemi == sr.het_or_hom_or_hemi)
            ).order_by(Sample.id.asc())]

            sr.sample_i = available_sample_ids.index(sr.id)  # the 1st sample will have sample_i = 0, etc.

            logging.info("%s-%s-%s-%s %s - computed sample_i: %s" % (
                sr.chrom, sr.pos, sr.ref, sr.alt, sr.het_or_hom_or_hemi, sr.sample_i))

        if sr.original_bam_path is None:
            sr.original_bam_path = lookup_original_bam_path(sr.sample_id)  # recompute the original bam path in case it's changed
            sr.save()

        if sr.sample_i >= MAX_SAMPLES_TO_SHOW_PER_VARIANT + BACKUP_SAMPLES_IN_CASE_OF_ERRORS:
            logging.info("%s-%s-%s-%s %s - sample_i too large. Skipping: %s" % (
                sr.chrom, sr.pos, sr.ref, sr.alt, sr.het_or_hom_or_hemi, sr.sample_i))
            sr.delete_instance()
            continue

        # make sure the variant is normalized
        if len(sr.ref) > 1 and len(sr.alt) > 1:
            #variant_records = list(
            #    Variant.select().where((Variant.chrom == sr.chrom) & (Variant.pos == sr.pos) & (Variant.ref == sr.ref) & (Variant.alt == sr.alt)))

            # normalize and update sample record and variant records
            (sr.chrom, sr.pos, sr.ref, sr.alt) = normalize(PYSAM_FASTA, str(sr.chrom), sr.pos, sr.ref, sr.alt)
            sr.variant_id = "%s-%s-%s-%s" % (sr.chrom, sr.pos, sr.ref, sr.alt)
            sr.save()
            #for v in variant_records:
            #    (v.chrom, v.pos, v.ref, v.alt, v.variant_id) = (sr.chrom, sr.pos, sr.ref, sr.alt, sr.variant_id)
            #    v.save()

        # check whether MAX_SAMPLES_TO_SHOW_PER_VARIANT have already been generated, and skip if yes
        if sr.sample_i >= MAX_SAMPLES_TO_SHOW_PER_VARIANT:
            # check whether enough samples have already been generated for this site
            num_samples_already_finished = len([s.id for s in Sample.select().where(
                (Sample.chrom == sr.chrom) & (Sample.pos == sr.pos) & (Sample.ref == sr.ref) & (Sample.alt == sr.alt) &
                (Sample.het_or_hom_or_hemi == sr.het_or_hom_or_hemi) & (Sample.started == 1 & Sample.finished == 1))])
            if num_samples_already_finished >= MAX_SAMPLES_TO_SHOW_PER_VARIANT:
                sr.started = 1
                sr.finished = 1
                sr.hc_error_code = 5000   
                sr.hc_error_text = "not needed"
                sr.save()
                continue


        try:
            run_haplotype_caller(
                sr.chrom,
                sr.pos,
                sr.ref,
                sr.alt,
                sr.het_or_hom_or_hemi,
                original_bam_path=sr.original_bam_path,
                original_gvcf_path=sr.original_gvcf_path,
                sample_id=sr.sample_id,
                sample_i=sr.sample_i,
                all_bam_output_dir=bam_output_dir,
                only_choose_samples=False,
            )
        except Exception as e:
            logging.error("%s-%s-%s-%s %s - error in run_haplotype_caller: %s" % (
                sr.chrom, sr.pos, sr.ref, sr.alt, sr.het_or_hom_or_hemi, e))
            traceback.print_exc()

    logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)
    logging.info("generate_HC_bams finished.") # at %s:%s" % (chrom, pos))


def create_sample_record_iterator(chrom=None, start_pos=None, end_pos=None):
    """Iterate over sample records that are marked as not-yet-finished in the sample table

    Args:
        chrom: chromosome
        start_pos: integer 1-based inclusive start position of genomic region
        end_pos: integer 1-based inclusive end position of genomic region
    Returns:
        Sample records
    """

    where_condition = (Sample.started == 0) & (Sample.finished == 0)
    where_condition = where_condition & (Sample.priority == 1)
    if chrom is not None:
        where_condition = where_condition & (Sample.chrom == chrom)
    if start_pos is not None:
        where_condition = where_condition & (Sample.pos >= start_pos)
    if end_pos is not None:
        where_condition = where_condition & (Sample.pos <= end_pos)

    while True:
        # get the next Sample to process
        randomized_sample_num = random.randint(1, 1000)
        unprocessed_samples = Sample.select().where(where_condition).order_by(
            Sample.sample_i.asc()).limit(randomized_sample_num)

        unprocessed_samples_list = list(unprocessed_samples)
        if len(unprocessed_samples_list) == 0:
            logging.info("Finished all samples. Exiting..")
            break

        logging.info("-----")
        #sql, params = unprocessed_samples.sql()
        #logging.info("query: %s \n rows retreived %s" % ( (sql % tuple(params)), len(unprocessed_samples_list)))

        current_sample = unprocessed_samples_list[-1]  # pick the last one
        logging.info("retrieving next sample id = %s: %s-%s-%s-%s %s" % (current_sample.id,
                                                                         current_sample.chrom,
                                                                         current_sample.pos,
                                                                         current_sample.ref,
                                                                         current_sample.alt,
                                                                         current_sample.het_or_hom_or_hemi))

        with db.atomic():
            query = Sample.update(started=1).where( (Sample.id == current_sample.id) & where_condition )
            rows_updated = query.execute()

        #sql, params = query.sql()
        #logging.info("query: %s \n rows updated: %s" % ( (sql % tuple(params)), rows_updated))

        if rows_updated == 0:
            sleep_interval = random.randint(1, 15)
            logging.info("%s-%s-%s-%s - sample claimed by another task. Skipping.. (sleep for %s sec)" % (
                current_sample.chrom, current_sample.pos, current_sample.ref, current_sample.alt, sleep_interval))
            time.sleep(sleep_interval) # sleep for a random time interval to avoid constant lock contension
            continue

        yield current_sample



if __name__ == "__main__":
    p = configargparse.getArgumentParser()
    p.add("--bam-output-dir", help="Where to output HC-reassembled bams", default=BAM_OUTPUT_DIR)
    p.add("--chrom", help="If specified, only process this chromosome", nargs="*")
    p.add("--start-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int)
    p.add("--end-pos", help="If specified, only process region in this interval (1-based inclusive coordinates)", type=int, default=10**10)

    p.add("--exit-after", metavar="MINUTES", help="This many minutes after starting, finish processing "
                                                  "the current sample and then exit", default=60*1, type=float)

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

    # db = init_db()  # commented out to avoid overloading database initially.
    db = _readviz_db

    try:
        from pyinstrument import Profiler
        profiling_enabled = True
    except:
        profiling_enabled = False

    profiling_enabled = False
    if profiling_enabled:
        profiler = Profiler() # or Profiler(use_signal=False), see below
        profiler.start()

    # process the samples
    logging.info("Processing remaining samples in sample table")
    if not args.chrom:
        chromosomes = [None]
    else:
        chromosomes = args.chrom

    for chrom in args.chrom:
        sample_record_iterator = create_sample_record_iterator(chrom=chrom, start_pos=args.start_pos, end_pos=args.end_pos)
        main(sample_iterator=sample_record_iterator, bam_output_dir=args.bam_output_dir, exit_after_minutes=args.exit_after)

    if profiling_enabled:
        profiler.stop()
        print(profiler.output_text(unicode=False, color=True))
