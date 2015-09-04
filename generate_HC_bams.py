"""
This script generates HC-reassembled bams. For each variant, it
selects which samples will be displayed and then runs HaplotypeCaller on those
samples. An unlimited number of instances of this script can be
run in parallel (for example in an array job) to speed up execution.
"""

# how many samples to show per het or hom-alt variant in the exac browser.
MAX_SAMPLES_TO_SHOW_PER_VARIANT = 5

# how many subdirectories to use to store the reassembled bams
NUM_OUTPUT_DIRECTORIES_L1 = 100
NUM_OUTPUT_DIRECTORIES_L2 = 10000

ACTIVE_REGION_PADDING = 300

import argparse
import collections
import datetime
import gzip
import os
import pysam
import re
import sqlite3
import subprocess
import sys
import tempfile
import time
import vcf

from postprocess_reassembled_bam import postprocess_bam
from minimal_representation import get_minimal_representation



MAX_ALLELE_LENGTH = 75

def compute_reassembled_bam_path(base_dir, chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom, sample_i, suffix=""):
    """Returns a reassembled bam output path"""
    output_subdir = "%02d/%04d" % (minrep_pos % NUM_OUTPUT_DIRECTORIES_L1,
                                   minrep_pos % NUM_OUTPUT_DIRECTORIES_L2)
    output_bam_filename = "chr%s-%s-%s-%s_%s%s%s.bam" % (
        chrom,
        minrep_pos,
        minrep_ref[:MAX_ALLELE_LENGTH],
        minrep_alt[:MAX_ALLELE_LENGTH],
        het_or_hom,
        sample_i,
        suffix)

    return os.path.join(base_dir, output_subdir, output_bam_filename)


def run(shell_cmd, verbose=False):
    """Utility method to print and run a shell command"""
    if verbose:
        print(shell_cmd)
    return os.system(shell_cmd)


def create_symlink(sample_bam_path, sample_gvcf_path, output_bam_path):
    output_dir = os.path.dirname(output_bam_path)
    if not os.path.isdir(output_dir):
        run("mkdir -p %s" % output_dir)

    symlink_path = output_bam_path.replace(".bam", "") + ".original.bam"
    if os.path.isfile(symlink_path):
        return
    run("ln -s %s %s" % (sample_bam_path, symlink_path))
    run("ln -s %s %s" % (sample_bam_path.replace(".bam", ".bai"), symlink_path + ".bai"))

    symlink_path = output_bam_path.replace(".bam", "") + ".original.gvcf"
    if os.path.isfile(symlink_path):
        return
    run("ln -s %s %s" % (sample_gvcf_path, symlink_path))
    run("ln -s %s %s" % (sample_gvcf_path + ".tbi", symlink_path + ".tbi"))


def main(args):
    """args: object returned by argparse.ArgumentParser.parse_args()"""

    # use exac info table to create in-memory mapping of sample_id => bam_path and sample_id => include
    sample_name_to_bam_path = {}
    sample_name_to_gvcf_path = {}   # original GVCF produced as part of the ExAC pipeline run
    sample_name_include_status = {}


    # db columns:  gvcf_call_mismatch, mismatch_type1, mismatch_type2, 
    tabix_file = pysam.TabixFile(filename=args.full_vcf, parser=pysam.asTuple())
    last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
    if args.chrom:
        vcf_iterator = tabix_file.fetch(args.chrom, args.start, args.end)  
    else:
        vcf_iterator = pysam.tabix_iterator(gzip.open(args.full_vcf), parser=pysam.asTuple())

    parse_vcf_row = create_vcf_row_parser(last_header_line, set(sample_name_include_status.keys()))

    # ok-to-be-public db - use this database to keep track of which variants have been completed
    db_table_columns = [
        "chrom text",
        "minrep_pos integer",
        "minrep_ref text",
        "minrep_alt text",

        "n_het integer",
        "n_hom_alt integer",
        "reassembled_bams_het text",  # comma-separated list of reassembled bam filenames for HET samples
        "reassembled_bams_hom text",  # comma-separated list of reassembled bam filenames for HOM samples

        "finished bool",

        #"start_date integer",
        #"finish_date integer",
        #"error text",
    ]


    temp_variants_db2_path = os.path.join("/tmp", "exac_v3_private_metadata%s.db" % db_filename_suffix) # use /tmp dir because SQLite doesn't work on NFS drives
    if os.path.isfile(temp_variants_db2_path):
        os.remove(temp_variants_db2_path)
    variants_db2 = sqlite3.connect(temp_variants_db2_path) #, isolation_level=False)
    variants_db2.execute("CREATE TABLE IF NOT EXISTS t(%s)" % ",".join(db2_table_columns))
    variants_db2.execute("CREATE UNIQUE INDEX IF NOT EXISTS variant_idx ON t(chrom, minrep_pos, minrep_ref, minrep_alt)")

    counters = collections.defaultdict(int)
    for row in vcf_iterator:
        counters["    sites"] += 1
        if counters["    sites"] % args.n_threads != args.thread_i:
            # extremely simple way to allow parallelization by running
            # multiple instances of this script: just skip sites
            # where site_counter % num_threads != thread_i
            continue

        chrom, pos, ref, alt_alleles, genotypes = parse_vcf_row(row)

        for alt_allele_index, alt in enumerate(alt_alleles):
            counters["   all_alleles"] += 1

            # minrep
            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(pos, ref, alt)

            print("############")
            print("SITE: %(chrom)s:%(minrep_pos)s-%(minrep_ref)s-%(minrep_alt)s" % locals())

            # skip variants that are already finished
            #(already_finished,) = variants_db.execute(
            #    "SELECT count(*) FROM t "
            #    "WHERE chrom=? AND minrep_pos=? AND minrep_ref=? AND minrep_alt=? and finished=1", (
            #    chrom, minrep_pos, minrep_ref, minrep_alt)).fetchone()
            #if already_finished:
            #    continue

            counters["  alleles_to_be_added_to_db"] += 1

            if args.use_calling_intervals:
                # look up the exac calling interval that overlaps this variant, so it can be passed to HaplotypeCaller (-L arg)
                regions = list(exac_calling_intervals_db.execute(
                    "SELECT chrom, start, end FROM intervals WHERE chrom=? AND start<=? AND end>=?", (chrom, minrep_pos, minrep_pos)))

                assert len(regions) != 0, "No region overlaps variant %s-%s: %s" % (chrom, minrep_pos, regions)
                assert len(regions) < 2, "Multiple regions overlap variant %s-%s: %s" % (chrom, minrep_pos, regions)

                region_chrom, region_start, region_end = regions[0]
                suffix = ""
            else:
                region_chrom, region_start, region_end = chrom, minrep_pos - ACTIVE_REGION_PADDING, minrep_pos + ACTIVE_REGION_PADDING
                suffix = "__+-%s" % ACTIVE_REGION_PADDING

            region = "%s:%s-%s" % (region_chrom, region_start, region_end)

            # print some stats
            if counters["   all_alleles"] % 10 == 0:
                if counters["   all_alleles"] % 10 == 0:
                    run("cp %s %s" % (temp_variants_db_path, args.bam_output_dir), verbose=True)
                    run("cp %s %s" % (temp_variants_db2_path, args.bam_output_dir), verbose=True)
                print("%s: %s.  %s  %s" % (
                    str(datetime.datetime.now()).split(".")[0],
                    ", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),
                    "-".join(map(str, (chrom, minrep_pos, minrep_ref, minrep_alt))),
                    region))

            # choose het and hom-alt samples to display
            chosen_samples = {}
            for het_or_hom in ["het", "hom"]:
                chosen_samples[het_or_hom] = choose_samples(
                    het_or_hom,
                    alt_allele_index + 1,  # add 1 since VCF genotypes count alt alleles starting from 1
                    genotypes,
                    sample_name_include_status,
                    sample_name_to_bam_path,
                    sample_name_to_gvcf_path)

            # skip variants where none of the samples were called het or hom-alt
            if len(chosen_samples["het"]) + len(chosen_samples["hom"]) == 0:
                print("No het or hom samples. Skipping this site...")
                counters[" hom_ref_alleles"] += 1
                continue

            # add variant to variants_db
            values = (chrom, minrep_pos, minrep_ref, minrep_alt,
                      len(chosen_samples["het"]),
                      len(chosen_samples["hom"]),
                      "", "", 0)
            question_marks = ",".join(["?"]*len(values))
            variants_db.execute("INSERT OR IGNORE INTO t VALUES ("+question_marks+")", values)
            variants_db.commit()

            # add variant to variants_db2
            values = (chrom, minrep_pos, minrep_ref, minrep_alt,
                      region_start,
                      region_end,
                      ",".join(map(lambda paths: paths[0], chosen_samples["het"])),
                      ",".join(map(lambda paths: paths[0], chosen_samples["hom"])),)
            question_marks = ",".join(["?"]*len(values))
            variants_db2.execute("INSERT OR IGNORE INTO t VALUES ("+question_marks+")", values)
            variants_db2.commit()

            counters["added_to_db"] += 1

            # run HaplotypeCaller to generated reassembled bam for each sample
            previously_seen_output_bam_paths = set()
            reassembled_bam_paths = collections.defaultdict(list)
            relative_reassembled_bam_paths = collections.defaultdict(list)
            for het_or_hom in ["het", "hom"]:
                chosen_samples_list = chosen_samples[het_or_hom]
                print("=======")
                print(str(len(chosen_samples_list)) + " chosen %(het_or_hom)s samples for %(chrom)s:%(minrep_pos)s-%(minrep_ref)s-%(minrep_alt)s: %(chosen_samples_list)s" % locals())
                for sample_i, (sample_bam_path, sample_gvcf_path) in enumerate(chosen_samples_list):
                    print("-----")
                    output_bam_path = compute_reassembled_bam_path(
                        args.bam_output_dir,
                        chrom,
                        minrep_pos,
                        minrep_ref,
                        minrep_alt,
                        het_or_hom,
                        sample_i,
                        suffix=suffix)

                    reassembled_bam_paths[het_or_hom].append(output_bam_path)

                    relative_output_bam_path = output_bam_path.replace(args.bam_output_dir+"/", "")
                    relative_reassembled_bam_paths[het_or_hom].append(relative_output_bam_path)

                    # sanity check - make sure bam path is not duplicate
                    if relative_output_bam_path in previously_seen_output_bam_paths:
                        print("ERROR: %s is not unique" % relative_output_bam_path)
                    else:
                        previously_seen_output_bam_paths.add(relative_output_bam_path)

                    if args.run_haplotype_caller:
                        if os.access(output_bam_path, os.R_OK):
                            print(("WARNING: reassembled bam already exists even though it's not marked "
                                "as finished in the database: %s. Will mark it as finished..") % output_bam_path)
                            continue

                        # run haplotype caller if reassembled bam doesn't exist yet
                        launch_haplotype_caller(
                            region,
                            args.fasta,
                            sample_bam_path,
                            output_bam_path)

                        new_gvcf_path = output_bam_path.replace(".bam", "")+".gvcf"
                        gvcfs_match = check_gvcf(sample_gvcf_path, new_gvcf_path, chrom, minrep_pos, minrep_ref, minrep_alt, het_or_hom)
                        if gvcfs_match:
                            if os.path.isfile(new_gvcf_path):
                                os.remove(new_gvcf_path)
                            if os.path.isfile(new_gvcf_path + ".idx"):
                                os.remove(new_gvcf_path + ".idx")

                        counters["bam_generated"] += 1

                    if args.create_links_to_original_bams or not gvcfs_match:
                        create_symlink(sample_bam_path, sample_gvcf_path, output_bam_path)

            if args.take_igv_screenshots:
                print("Taking igv screenshots")
                import igv_api
                #igv_jar_path="/home/unix/mlek/bin/IGV_plotter/IGV_2.0.35/igv.jar"
                igv_jar_path="~/bin/igv-2.3.57/igv.jar"
                r = igv_api.IGVCommandLineRobot(verbose=True,
                                                igv_window_width=1200,
                                                igv_window_height=1600,
                                                igv_jar_path=igv_jar_path)

                for het_or_hom in ["het", "hom"]:
                    r.new_session()
                    r.max_panel_height(2000)
                    original_bams = map(lambda paths: paths[0], chosen_samples[het_or_hom])
                    for original_bam, reassembled_bam in zip(original_bams, reassembled_bam_paths[het_or_hom]):
                        print("%s vs %s" % (original_bam, reassembled_bam))
                    r.load(list(sum(zip(original_bams, reassembled_bam_paths[het_or_hom]), ())))
                    r.load(["./gencode.v19.annotation.gtf.gz",
                            "./exome_calling_regions.v1.bed.gz"])
                    r.goto("%s:%s-%s" % (region_chrom, minrep_pos - 205, minrep_pos + 195))
                    output_png_path = os.path.join(os.path.dirname(output_bam_path), os.path.basename(output_bam_path).split("_")[0]+"_"+het_or_hom+suffix+".png")
                    r.screenshot(output_png_path)
                    r.exit_igv()
                    r.execute()

            # update variants_db
            variants_db.execute("UPDATE t SET "
                "finished=1, reassembled_bams_het=?, reassembled_bams_hom=?"
                "WHERE "
                "chrom=? AND minrep_pos=? AND minrep_ref=? AND minrep_alt=?", (
                    ",".join(relative_reassembled_bam_paths["het"]),
                    ",".join(relative_reassembled_bam_paths["hom"]),
                    chrom, minrep_pos, minrep_ref, minrep_alt))
            variants_db.commit()

    variants_db.close()

    print("Finished.")
    run("cp %s %s" % (temp_variants_db_path, args.bam_output_dir), verbose=True)
    run("cp %s %s" % (temp_variants_db2_path, args.bam_output_dir), verbose=True)



if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--info-table", help="Path of ExAC info table",
        default='/humgen/atgu1/fs03/lek/resources/ExAC/ExAC.r0.3_meta_Final.tsv')
    p.add_argument("--full-vcf", help="Path of the ExAC full vcf with genotypes",
        default='/humgen/atgu1/fs03/konradk/exac/gqt/exac_all.vcf.gz')
    p.add_argument("-R", "--fasta", help="Reference genome hg19 fasta",
        default="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
    p.add_argument("--bam-output-dir", help="Where to output HC-reassembled bams",
        default='/broad/hptmp/exac_readviz_backend/')
    p.add_argument("--exac-calling-intervals", help="ExAC calling regions .intervals file",
        default="/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list")
    p.add_argument("-i", "--thread-i", help="Thread number (must be between 1 and the value of --n-threads)", default=1, type=int)
    p.add_argument("-n", "--n-threads", help="Total number of threads", default=1, type=int)
    p.add_argument("--take-igv-screenshots", help="Whether to take IGV snapshots", action="store_true")
    p.add_argument("--run-haplotype-caller", help="Whether to run haplotype caller", action="store_true")
    p.add_argument("--create-links-to-original-bams", help="Whether to create symlinks to the original bams for each sample for debugging", action="store_true")
    #p.add_argument("--use-calling-intervals", help="Whether to pass the ExAC calling region to -L when running HaplotypeCaller", action="store_true")
    p.add_argument("chrom", nargs="?", help="If specified, only data for this chromosome will be added to the database.")
    args = p.parse_args()

    args.use_calling_intervals = True # use them by default

    # print out settings
    print("Running with settings: ")
    for argname in filter(lambda n: not n.startswith("_"), dir(args)):
        print("   %s = %s" % (argname, getattr(args, argname)))
    print("\n")

    # validate args
    i = 0
    while True:
        try:
            assert os.path.isfile(args.info_table), "Couldn't find: %s" % args.info_table
            assert os.path.isfile(args.full_vcf), "Couldn't find: %s" % args.full_vcf
            assert os.path.isfile(args.full_vcf+".tbi"), "Couldn't find: %s.tbi" % args.full_vcf
            assert os.path.isdir(args.bam_output_dir), "Couldn't find: %s" % args.bam_output_dir
            assert os.path.isfile(args.exac_calling_intervals), "Couldn't find: %s" % args.exac_calling_intervals

            assert args.thread_i > 0, "Invalid --thread_i arg (%s) or --n-threads arg (%s)" % (args.thread_i, args.n_threads)
            assert args.thread_i <= args.n_threads, "Invalid --thread_i arg (%s) or --n-threads arg (%s)" % (args.thread_i, args.n_threads)
            break
        except AssertionError:
            if i >= 5:
                raise
            i += 1  # retry logic
            time.sleep(1)


    if args.chrom:
        args.chrom = args.chrom.replace("chr", "").upper()
        if args.chrom not in list(map(str, range(1, 23))) + ["X", "Y", "M", "MT"]:
            p.error("Invalid chromosome: " + args.chrom)
        print("Thread #%s out of %s will run on chr%s" % (args.thread_i, args.n_threads, args.chrom))
    else:
        print("Thread #%s out of %s will run on all chromosomes" % (args.thread_i, args.n_threads))

    args.thread_i = args.thread_i - 1  # convert to a 0-based count
    main(args)
