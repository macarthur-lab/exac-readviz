"""
This module parses the calling intervals file into a database table (if
this hasn't been done already). The table is indexed by (chrom, start, end).
"""
import logging
import os
import peewee
import sys

# utility functions
def parse_exac_calling_intervals(exac_calling_intervals_path):
    """Parses the exac calling regions .intervals file into a database
    'ExacCallingInterval' table that has columns:
        chrom, start_pos, end_pos, strand, target_name.

    Args:
      exac_calling_intervals_path: A path like "/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list"
    """
    logging.info("Loading exac calling interals from %s" % exac_calling_intervals_path)

    # check if already done
    with open(exac_calling_intervals_path) as f:
        lines_in_file = sum(1 for line in f if not line.startswith("@"))
    records_in_db = ExacCallingInterval.select().count()
    if records_in_db > 0:
        message = "%(records_in_db)d out of %(lines_in_file)d intervals from %(exac_calling_intervals_path)s are already in the database" % locals()
        if lines_in_file == records_in_db:
            logging.info("All " + message + ". Exiting...")
            return
        else:
            logging.info(message)

    # parse exac_calling_intervals_file into the database table
    with open(exac_calling_intervals_path) as f:
        for line in f:
            # skip header
            if line.startswith("@"):
                continue

            # parse fields: chrom, start_pos, end_pos, strand, target_name
            fields = line.strip("\n").split("\t")
            ExacCallingInterval.get_or_create(
                chrom=fields[0].replace("chr", ""),
                start=fields[1],
                end=fields[2],
                strand=fields[3],
                name=fields[4])

    # print some stats
    for row in ExacCallingInterval.raw("select count(*) as n, min(end-start) as min_size, max(end-start) as max_size, avg(end-start) as avg_size from " + ExacCallingInterval.__name__.lower()):
        logging.info("Loaded %s intervals (range: %s to %s, mean: %s) into %s" % (row.n, row.min_size, row.max_size, int(row.avg_size), calling_intervals_db))
    logging.info("Done loading")


def get_overlapping_calling_interval(chrom, pos):
    """Returns the calling interval that contains the given chrom, pos
    (1-based, inclusive)
    """
    intervals = list(
        ExacCallingInterval.select("chrom", "start", "end").where(
            chrom=chrom, start <= pos, end >= pos))

    if not intervals:
        raise ValueError("No region overlaps variant %s-%s" % (chrom, pos))


    assert len(intervals) < 2, "Multiple regions overlap variant %s-%s: %s" % (
        chrom, pos, intervals)

    i = intervals[0]
    return i.chrom, i.start, i.end



calling_intervals_db = peewee.MySQLDatabase('exac_readviz', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307)
calling_intervals_db.connect()

# define database model for storing the calling intervals
class ExacCallingInterval(peewee.Model):
    chrom = peewee.CharField(max_length=5, null=False, index=True)  # without 'chr'
    start = peewee.IntegerField()
    end = peewee.IntegerField()
    strand = peewee.CharField(max_length=1)
    name = peewee.CharField(max_length=500)

    class Meta:
        database = calling_intervals_db

indexes = (
    (('chrom', 'start', 'end'), True),  # True means unique index
)


import configargparse

from utils.database import create_table

p = configargparse.getArgumentParser()
p.add("--exac-calling-intervals", help="ExAC calling regions .intervals file",
      default="/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list")
args = p.parse_args()

# create the database table and indexes if they don't exist yet
if not ExacCallingInterval.table_exists():
    # parse the exac calling intervals file into a database table
    assert os.path.isfile(args.exac_calling_intervals), \
        "Couldn't find: %s" % args.exac_calling_intervals

    create_table(calling_intervals_db, ExacCallingInterval, indexes, safe=True)

    parse_exac_calling_intervals(args.exac_calling_intervals)
