
# inputs: directories
# read in each bam, check that it's
# output: combined, indexed bam
# for each input directory

import argparse
import glob
import logging
import os
import peewee
from playhouse import shortcuts
import pysam
import traceback

from utils.database import Variant
from utils.constants import BAM_OUTPUT_DIR, PICARD_JAR_PATH


def run(cmd):
    logging.info(cmd)
    os.system(cmd)


def bam_path_to_fields(bam_path):
    # for example: /read_viz/22/5822/chr22-46615822-A-G_het0.bam
    return os.path.basename(bam_path).replace(".bam", "").replace('_', '-').split('-')


def bam_path_to_dict(bam_path):
    # for example: /read_viz/22/5822/chr22-46615822-A-G_het0.bam
    return dict(zip(['chrom', 'pos', 'ref', 'alt', 'het_or_hom'], bam_path_to_fields(bam_path)))


def bam_path_to_read_group_id(bam_path):
    return os.path.basename(bam_path.replace("chr", "").replace(".bam", ""))


def combine_bams(output_dir, chrom, position_hash, force=False):
    """Generates the combined bam.
    Args:
        force: proceed even if the .bam and .db are already there on disk.
    """

    hash_dir = "%03d" % (position_hash % 1000)

    obam_path = os.path.join(output_dir, "combined_bams", chrom, "_tmp.combined_chr%s_%s.bam" % (chrom, hash_dir))
    sorted_bam_path = obam_path.replace("_tmp.", "")

    sqlite_db_path = os.path.basename(obam_path.replace(".bam", ".db"))
    final_sqlite_db_path = sorted_bam_path.replace(".bam", ".db")

    # check if combine_bam output files already exist on disk, and skip if yes
    if not force and os.path.isfile(sorted_bam_path) and os.path.isfile(final_sqlite_db_path):
        sorted_bam_file_size = os.path.getsize(sorted_bam_path)
        if sorted_bam_file_size > 1000:
            logging.info("%s found on disk. size=%s. Skipping..." % (sorted_bam_path, sorted_bam_file_size))
            return

    expected_variants = [v for v in peewee.RawQuery(Variant,  (
        "select * from variant where chrom='%s' and pos %%%% 1000 = %s and finished=1") % (
        chrom, position_hash)).execute()]

    num_expected_to_be_available_bams = sum([v.n_available_samples for v in expected_variants if v.n_available_samples])
    logging.info("num_expected_to_be_available_bams = %(num_expected_to_be_available_bams)s (from the database)" % locals())

    # check number of expected variants
    expected_to_be_available_bam_paths = set()
    for v in expected_variants:
        if v.n_available_samples:
            expected_to_be_available_bam_paths.update(v.readviz_bam_paths.split("|"))

    # check how many are actually on disk
    actually_available_bam_paths = set(glob.glob(os.path.join(output_dir, chrom, hash_dir, 'chr*.bam')))

    num_actually_available_bams = len(actually_available_bam_paths)
    logging.info("num_actually_available_bams = %(num_actually_available_bams)s (from looking on disk)" % locals())

    if num_expected_to_be_available_bams > len(actually_available_bam_paths):
        logging.error("ERROR: %(num_expected_to_be_available_bams)s > %(num_actually_available_bams)s " % locals())

    # make sure that all of the expected_to_be_available_bam_paths are actually available
    actually_available_bam_paths = set(map(lambda p: p.replace(output_dir, ""), actually_available_bam_paths))

    if len(expected_to_be_available_bam_paths - actually_available_bam_paths) > 0:
        message = "ERROR: expected_to_be_available_bam_paths - actually_available_bam_paths: %s" % str(
            expected_to_be_available_bam_paths - actually_available_bam_paths)
        logging.error(message)
        raise Exception(message)

    logging.info("Processing %d expected (%d existing) bams in %s/%s/%s" % (
        len(expected_variants), len(actually_available_bam_paths), output_dir, chrom, hash_dir))

    # actually combine the bams. As confirmed above, the intersection of the 2 sets is identical to the expected_to_be_available_bam_paths set.
    print("Expected bams: " + str(list(expected_to_be_available_bam_paths)[0:5]))
    print("Available bams: " + str(list(actually_available_bam_paths)[0:5]))
    ibam_paths = list(expected_to_be_available_bam_paths & actually_available_bam_paths)

    logging.info("Combining %s bams into %s" % (len(ibam_paths), obam_path))
    if len(ibam_paths) > 0:
        # sort bams by position so that the reads in the combined file are roughly in sorted order
        sorted_bam_paths = sorted(ibam_paths, key=lambda ibam_path: int(bam_path_to_dict(ibam_path)['pos']))

        read_group_ids = map(bam_path_to_read_group_id, sorted_bam_paths)
        read_groups = [{'ID': read_group_id, "SM": 0} for read_group_id in read_group_ids]

        #ibam_paths = []
        obam = None
        for ibam_path in ibam_paths:
            try:
                ibam = pysam.AlignmentFile(os.path.join(output_dir, ibam_path), "rb")

                if obam is None:
                    header = {
                        'HD': {'VN': '1.4'},  #, 'SO': 'coordinate'},
                        'SQ': ibam.header['SQ'],
                        'RG': read_groups,
                    }

                    obam = pysam.AlignmentFile(obam_path, "wb", header=header)

                # iterate over the reads
                rg_tag = (('RG', bam_path_to_read_group_id(ibam_path)), )
                for r in ibam:
                    r.tags = rg_tag
                    obam.write(r)

                ibam.close()

            except (IOError, ValueError) as e:
                logging.error("ERROR: %s", e)
                logging.error(traceback.format_exc())
        if obam is not None:
            obam.close()


        # sort the file
        logging.info("Running picard SortSam:")
        picard_jar = PICARD_JAR_PATH

        run(("java -jar %(picard_jar)s SortSam VALIDATION_STRINGENCY=LENIENT "
            "I=%(obam_path)s O=%(sorted_bam_path)s SO=coordinate CREATE_INDEX=true") % locals())
        run("rm %(obam_path)s" % locals())

    # create a combined_chr<chrom>_<hash>.db sqlite database that the website code can use to understand whats in the combined file
    logging.info("Populating sqlite database: " + sqlite_db_path)
    if os.path.isfile(sqlite_db_path):
        run("rm " + sqlite_db_path)

    sqlite_db = peewee.SqliteDatabase(sqlite_db_path, autocommit=False)
    class t(Variant):
        class Meta:
            database = sqlite_db
            indexes = (
                (('chrom', 'pos', 'ref', 'alt', 'het_or_hom'), True), # True means unique index
            )

    t.create_table(fail_silently=True)

    # copy the records from the Variant table used by generate_HC_bams.py
    sqlite_db.connect()
    with sqlite_db.atomic():
        for v in expected_variants:  #Variant.select().where(Variant.finished==1).dicts():
            d = shortcuts.model_to_dict(v)

            # delete readviz_bam_paths as they're no longer relevant because the data from these is being combined into one bam file
            del d['readviz_bam_paths']

            t.insert(**d).execute()
    sqlite_db.close()

    run("mv %s %s" % (sqlite_db_path, final_sqlite_db_path))
    logging.info("-- interval finished --")  # detected by parallelize.py to mark this as done


if __name__ == "__main__":
    p = argparse.ArgumentParser("Generates combined bams")
    p.add_argument("-d", "--output-dir", help="the top-level output directory", default=BAM_OUTPUT_DIR)
    p.add_argument("-f", "--force", help="Regenerate combined .bam and sqlite .db even they already exist", action="store_true")
    p.add_argument("--chrom", help="optional chromosome", required=True)
    g = p.add_argument_group()
    g.add_argument("-k", "--position-hash",
       help="bams are divided between directories with names 000 through 999. "
            "This should be a number between 0 and 999 which specifies which of "
            "these directories to process.", type=int)
    g.add_argument("-k1", "--start-pos",
       help="bams are divided between directories with names 000 through 999. "
            "This should be a number between 0 and 999 which specifies the start of "
            "a range of these directories to process.", type=int)
    g.add_argument("-k2", "--end-pos", type=int,
       help="bams are divided between directories with names 000 through 999. "
            "This should be a number between 0 and 999 which specifies the end of "
            "a range of these directories to process.")
    args = p.parse_args()

    if args.position_hash is not None:
        combine_bams(args.output_dir, args.chrom, args.position_hash, force=args.force)
    elif args.start_pos is not None and args.end_pos is not None:
        for position_hash in range(args.start_pos, args.end_pos+1):
            logging.info("-------")
            combine_bams(args.output_dir, args.chrom, position_hash, force=args.force)
    else:
        p.error("Must specify -k or both -k1 and -k2")

#CREATE TABLE t(
# chrom text,
# minrep_pos integer,
# minrep_ref text,
# minrep_alt text,
# n_het integer,
# n_hom_alt integer,
# reassembled_bams_het text,
# reassembled_bams_hom text,
# finished bool);

# CREATE UNIQUE INDEX variant_idx ON t(
#   chrom,
#   minrep_pos,
#   minrep_ref,
#   minrep_alt);
