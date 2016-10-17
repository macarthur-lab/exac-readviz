
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

from utils.haplotype_caller import compute_output_bam_path
from utils.database import Variant, Sample, _SharedVariantPositionFields
from utils.constants import BAM_OUTPUT_DIR, PICARD_JAR_PATH, MAX_SAMPLES_TO_SHOW_PER_VARIANT


def run(cmd):
    logging.info(cmd)
    os.system(cmd)


def bam_path_to_fields(bam_path):
    # for example: /read_viz/22/5822/chr22-46615822-A-G_het0.bam
    return os.path.basename(bam_path).replace(".bam", "").replace('_', '-').split('-')


def bam_path_to_dict(bam_path):
    # for example: /read_viz/22/5822/chr22-46615822-A-G_het0.bam
    return dict(zip(['chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi'], bam_path_to_fields(bam_path)))


def bam_path_to_read_group_id(bam_path):
    return os.path.basename(bam_path.replace("chr", "").replace(".bam", ""))


def choose_samples_to_combine(variants_to_process):
    all_chosen_bam_paths_in_db = set()
    for v in variants_to_process:
        successful_samples = list(Sample.select().where(
                (Sample.chrom == v.chrom) &
                (Sample.pos == v.pos) &
                (Sample.ref == v.ref) &
                (Sample.alt == v.alt) &
                (Sample.het_or_hom_or_hemi == v.het_or_hom_or_hemi) &
                (Sample.hc_succeeded == 1) &
                (Sample.started == 1) &
                (Sample.finished == 1)
        ).order_by(Sample.id.asc()))

        # choose the 1st min(n_expected_samples, MAX_SAMPLES_TO_SHOW_PER_VARIANT) samples
        num_samples_to_chose = min(v.n_expected_samples, MAX_SAMPLES_TO_SHOW_PER_VARIANT)

        chosen_samples = successful_samples[0:num_samples_to_chose]

        for i, s in enumerate(chosen_samples):
            if s.sample_i == i:
                s.output_bam_path2 = s.output_bam_path
            else:
                s.output_bam_path2 = compute_output_bam_path(
                        s.chrom, s.pos, s.ref, s.alt, s.het_or_hom_or_hemi, i)
                print(s.output_bam_path + " => NEW OUTPUT PATH: " + s.output_bam_path2)

        readviz_bam_paths_for_variant = [s.output_bam_path2 for s in chosen_samples]

        if len(chosen_samples) < num_samples_to_chose:
            logging.error("%s-%s-%s-%s %s - ERROR: expected %s samples. Found only %s successful samples in database: %s" % (
                v.chrom, v.pos, v.ref, v.alt, v.het_or_hom_or_hemi,
                num_samples_to_chose, len(chosen_samples), ", ".join(readviz_bam_paths_for_variant)))

        if len(readviz_bam_paths_for_variant) > len(set(readviz_bam_paths_for_variant)):
            logging.error("%s-%s-%s-%s %s - ERROR: %s duplicate readviz bam paths: %s" % (
                v.chrom, v.pos, v.ref, v.alt, v.het_or_hom_or_hemi,
                len(readviz_bam_paths_for_variant) - len(set(readviz_bam_paths_for_variant)),
                str(readviz_bam_paths_for_variant)))

        # update variant table
        v.n_available_samples = len(chosen_samples)
        v.readviz_bam_paths = "|".join(readviz_bam_paths_for_variant)
        v.save()

        all_chosen_bam_paths_in_db.update(readviz_bam_paths_for_variant)

    return all_chosen_bam_paths_in_db


def generate_combined_bam(base_dir, reassembled_bam_paths, temp_combined_bam_path, combined_bam_path):
    logging.info("combining %s bams into %s" % (len(reassembled_bam_paths), combined_bam_path))

    # sort bams by position so that the reads in the combined file are roughly in sorted order
    sorted_bam_paths = sorted(reassembled_bam_paths, key=lambda ibam_path: int(bam_path_to_dict(ibam_path)['pos']))

    read_group_ids = map(bam_path_to_read_group_id, sorted_bam_paths)
    read_groups = [{'ID': read_group_id, "SM": 0} for read_group_id in read_group_ids]

    #ibam_paths = []
    obam = None
    for reassembled_bam_path in reassembled_bam_paths:
        # try reading in the reassembled bam and adding it to the combined bam
        try:
            ibam = pysam.AlignmentFile(os.path.join(base_dir, reassembled_bam_path), "rb")

            if obam is None:
                header = {
                    'HD': {'VN': '1.4'},  #, 'SO': 'coordinate'},
                    'SQ': ibam.header['SQ'],
                    'RG': read_groups,
                }

                obam = pysam.AlignmentFile(temp_combined_bam_path, "wb", header=header)

            # iterate over the reads
            rg_tag = (('RG', bam_path_to_read_group_id(reassembled_bam_path)), )
            for r in ibam:
                r.tags = rg_tag
                obam.write(r)

            ibam.close()

        except (IOError, ValueError) as e:
            logging.error("ERROR on file %s: %s", reassembled_bam_path, e)
            logging.error(traceback.format_exc())
    if obam is not None:
        obam.close()


    # sort the file
    logging.info("Running picard SortSam:")
    picard_jar = PICARD_JAR_PATH

    run(("java -jar %(picard_jar)s SortSam VALIDATION_STRINGENCY=LENIENT "
         "I=%(temp_combined_bam_path)s O=%(combined_bam_path)s SO=coordinate CREATE_INDEX=true") % locals())
    run("rm %(temp_combined_bam_path)s" % locals())
    
    bai_path = combined_bam_path.replace(".bam", ".bai")
    run("cp %(bai_path)s %(combined_bam_path)s.bai" % locals())  # copy the .bai file to .bam.bai since this is what IGV.js looks for


def generate_sqlite_db(variants_to_process, temp_sqlite_db_path, sqlite_db_path):
    logging.info("populating sqlite database: " + temp_sqlite_db_path)
    if os.path.isfile(temp_sqlite_db_path):
        run("rm -f " + temp_sqlite_db_path)

    sqlite_db = peewee.SqliteDatabase(temp_sqlite_db_path, autocommit=False)
    class t(_SharedVariantPositionFields):
        n_expected_samples = peewee.IntegerField(index=True, null=True)
        n_available_samples = peewee.IntegerField(index=True, null=True)

        class Meta:
            database = sqlite_db
            indexes = (
                (('chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi'), True), # True means unique index
            )

    t.create_table(fail_silently=True)

    # copy the records from the Variant table used by generate_HC_bams.py
    sqlite_db.connect()
    with sqlite_db.atomic():
        for v in variants_to_process:  #Variant.select().where(Variant.finished==1).dicts():
            #shortcuts.model_to_dict(v)

            d = {
                'chrom': v.chrom, 'pos': v.pos, 'ref': v.ref, 'alt': v.alt, 'het_or_hom_or_hemi': v.het_or_hom_or_hemi,
                'n_expected_samples': v.n_expected_samples,
                'n_available_samples': v.n_available_samples,
            }

            # delete readviz_bam_paths as they're no longer relevant because the data from these is being combined into one bam file
            print("INSERTING " + str(d))
            t.insert(**d).execute()
    sqlite_db.close()

    run("mv %s %s" % (temp_sqlite_db_path, sqlite_db_path))


def combine_bams(output_dir, temp_dir, chrom, position_hash, force=False):
    """Generates the combined bam.
    Args:
        output_dir: the top level directory where files are stored.
        temp_dir: a local non-NFS directory that can be used for creating /
            modifying a SQLite database. This avoids SQLite incompatibility with
            network drives.
        chrom: chromosome for which to combine bams
        position_hash: bams are divided between directories with names 000
            through 999. This should be a number between 0 and 999 which
            specifies which of these directories to process.
        force: proceed even if the .bam and .db are already there on disk.
    """

    # create combined_bams output sub-directory
    hash_dir = "%03d" % (position_hash % 1000)

    temp_combined_bam_path = os.path.join(output_dir, "combined_bams", chrom, "_tmp.combined_chr%s_%s.bam" % (chrom, hash_dir))
    combined_bam_path = temp_combined_bam_path.replace("_tmp.", "")

    if not os.path.isdir(os.path.dirname(temp_combined_bam_path)):
        run("mkdir -m 777 -p " + os.path.dirname(temp_combined_bam_path))

    temp_sqlite_db_path = os.path.join(temp_dir, os.path.basename(temp_combined_bam_path.replace(".bam", ".db")))
    sqlite_db_path = combined_bam_path.replace(".bam", ".db")

    # create iterator over variants in this bin
    variants_to_process = [v for v in peewee.RawQuery(Variant,  (
        "select * from variant where chrom='%s' and pos %%%% 1000 = %s"
    ) % (chrom, position_hash)).execute()]

    # choose the samples to combine, and get their reassembled bam paths
    all_chosen_bam_paths_in_db = choose_samples_to_combine(variants_to_process)

    # check if combine_bam output files already exist on disk, and skip if yes
    if not force and os.path.isfile(combined_bam_path) and os.path.isfile(sqlite_db_path):
        try:
            ibam = pysam.AlignmentFile(combined_bam_path, "rb")
            num_read_groups = len(ibam.header['RG'])
            ibam.close()
            if num_read_groups == len(all_chosen_bam_paths_in_db):
                logging.info("%s found on disk. size=%s read groups. Skipping..." % (combined_bam_path, num_read_groups))
                return
        except (IOError, ValueError) as e:
            logging.warning("WARNING: couldn't read combined file %s: %s", combined_bam_path, e)
            logging.warning(traceback.format_exc())
            logging.warning("Will regenerate it..")

    # check that all_chosen_bam_paths_in_db exist on disk
    all_available_bam_paths_on_disk = set(glob.glob(os.path.join(output_dir, chrom, hash_dir, 'chr*.bam')))
    all_available_bam_paths_on_disk = set(map(lambda p: p.replace(output_dir, '').strip('/'), list(all_available_bam_paths_on_disk)))

    all_chosen_bam_paths_on_disk = all_available_bam_paths_on_disk & all_chosen_bam_paths_in_db

    if len(all_chosen_bam_paths_on_disk) < len(all_chosen_bam_paths_in_db):
        logging.info("ERROR: found only %s out of %s reassembled bams on disk in %s/%s/%s" % (
            len(all_chosen_bam_paths_on_disk), len(all_chosen_bam_paths_in_db), output_dir, chrom, hash_dir))
    else:
        logging.info("all %s out of %s reassembled bams found on disk in %s/%s/%s" % (
            len(all_chosen_bam_paths_on_disk), len(all_chosen_bam_paths_in_db), output_dir, chrom, hash_dir))

    if len(all_chosen_bam_paths_on_disk) > 0:
        generate_combined_bam(base_dir=output_dir, reassembled_bam_paths=list(map(str, all_chosen_bam_paths_on_disk)),
                              temp_combined_bam_path=temp_combined_bam_path, combined_bam_path=combined_bam_path)

    # create a combined_chr<chrom>_<hash>.db sqlite db where the website code can look up the number of expected and
    # available readviz tracks for each variant in this bin
    generate_sqlite_db(variants_to_process, temp_sqlite_db_path, sqlite_db_path)

    logging.info("-- interval finished --")  # detected by parallelize.py to mark this as done


if __name__ == "__main__":
    p = argparse.ArgumentParser("Generates combined bams")
    p.add_argument("--chrom", help="chromosome", required=True)
    p.add_argument("-d", "--output-dir", help="the top-level output directory", default=BAM_OUTPUT_DIR)
    p.add_argument("-t", "--non-nfs-temp-dir", help="local non-NFS-mounted temp directory to use for sqlite oprations", default="/tmp")
    p.add_argument("-f", "--force", help="regenerate combined .bam and sqlite .db even they already exist", action="store_true")
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
        combine_bams(args.output_dir, args.non_nfs_temp_dir, args.chrom, args.position_hash, force=args.force)
    elif args.start_pos is not None and args.end_pos is not None:
        for position_hash in range(args.start_pos, args.end_pos+1):
            logging.info("-------")
            combine_bams(args.output_dir, args.non_nfs_temp_dir, args.chrom, position_hash, force=args.force)
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
