"""
This module parses the calling intervals file into a database table (if
this hasn't been done already). The table is indexed by (chrom, start, end).
"""
import itertools
import logging
import os
import peewee
import sys

calling_intervals_db = peewee.MySQLDatabase('exac_readviz', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307)
#calling_intervals_db = peewee.MySQLDatabase('exac_readviz2', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307, autocommit=False)
#calling_intervals_db = peewee.SqliteDatabase('exac_readviz.db', autocommit=False)
#calling_intervals_db = peewee.SqliteDatabase(':memory:')
calling_intervals_db.connect()

#logging.getLogger('peewee').setLevel(logging.DEBUG)

# utility functions
def parse_exac_calling_intervals(exac_calling_intervals_path):
    """Parses the exac calling regions .intervals file into a database
    'ExacCallingInterval' table that has columns:
        chrom, start_pos, end_pos, strand, target_name.

    Args:
      exac_calling_intervals_path: A path like "/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list"
    """

    # check if already done
    with open(exac_calling_intervals_path) as f:
        lines_in_file = sum(1 for line in f if not line.startswith("@"))
    records_in_db = ExacCallingInterval.select().count()
    if records_in_db > 0:
        message = "%(records_in_db)d out of %(lines_in_file)d intervals were previously loaded from %(exac_calling_intervals_path)s" % locals()
        if lines_in_file == records_in_db:
            logging.info("All " + message)
            return
        else:
            logging.info("Only " + message)
    else:
        logging.info("Loading exac calling intervals from %s" % exac_calling_intervals_path)

    # parse exac_calling_intervals_file into the database table
    with open(exac_calling_intervals_path) as f:
        previous_chrom = None
        previous_start_pos = None
        line_counter = 0

        #calling_intervals_db.begin()
        for line in f:
            # skip header
            if line.startswith("@"):
                continue

            line_counter += 1
            if line_counter % 10000 == 0:
                logging.info("%d lines processed" % line_counter)
                #calling_intervals_db.commit()

            # parse fields: chrom, start_pos, end_pos, strand, target_name
            fields = line.strip("\n").split("\t")
            chrom = fields[0].replace("chr", "")
            start_pos, end_pos = map(int, fields[1:3])
            strand, name = fields[3:5]

            if previous_chrom != chrom:
                # reset
                previous_start_pos = start_pos
                previous_chrom = chrom

            # ensure intervals are in order by start-pos so that table's id can
            # be used to retrieve next and previous intervals
            assert previous_start_pos <= start_pos, \
                "Intervals in %s are out of order: %s is before %s" % (
                    exac_calling_intervals_path,
                    (previous_chrom, previous_start_pos),
                    (chrom, start_pos),
                )

            ExacCallingInterval.get_or_create(
                chrom=chrom,
                start=start_pos,
                end=end_pos,
                strand=strand,
                name=name)
    #calling_intervals_db.commit()

    # print some stats
    for row in ExacCallingInterval.raw("select count(*) as n, min(end-start) as min_size, max(end-start) as max_size, avg(end-start) as avg_size from " + ExacCallingInterval.__name__.lower()):
        logging.info("Loaded %s intervals (range: %s to %s, mean: %s) into %s" % (row.n, row.min_size, row.max_size, int(row.avg_size), calling_intervals_db))
    logging.info("Done loading")


def _get_interval(chrom, pos):
    """Utility method - returns the ExacCallingInterval object that overlaps the
    given chrom, pos (1-based, inclusive)
    """
    intervals = list(
        ExacCallingInterval.select().where(
            (ExacCallingInterval.chrom == chrom) & (ExacCallingInterval.start <= pos) & (ExacCallingInterval.end >= pos)))

    if not intervals:
        raise ValueError("Variant %s:%s is not overlapped by any calling intervals" % (chrom, pos))

    assert len(intervals) < 2, "Multiple calling intervals overlap variant %s:%s - %s" % (
        chrom, pos, intervals)

    return intervals[0]


def get_overlapping_calling_interval(chrom, pos):
    """Returns the (chrom, start end) for the calling interval that overlaps
    the given chrom, pos (1-based, inclusive).
    """
    i = _get_interval(chrom, pos)
    return i.chrom, i.start, i.end


def get_adjacent_calling_intervals(chrom, pos, n_left=1, n_right=1):
    """Returns the calling interval that overlaps the given chrom: pos (1-based,
    inclusive), and also up to n_left intervals to the left of this interval
    (in order by chrom & start position), as well as up to n_right intervals
    to the right of the interval. If the overlapping interval is at the very
    beginning or very end of a chromosome, fewer than n_left or n_right intervals
    may be included.

    Args:
      chrom: variant chrom
      pos: variant pos
      n_left: returned list will include this many intervals to the left of the
        interval that overlaps the given variant
      n_right: returned list will include this many intervals to the right of the
        interval that overlaps the given variant

    Return:
       list of intervals represented as (chrom, start, stop) tuples
    """
    i = _get_interval(chrom, pos)
    left_intervals = ExacCallingInterval.select().where(
        ExacCallingInterval.chrom == chrom,
        ExacCallingInterval.id >= i.id - n_left,
        ExacCallingInterval.id < i.id)
    right_intervals = ExacCallingInterval.select().where(
        ExacCallingInterval.chrom == chrom,
        ExacCallingInterval.id > i.id,
        ExacCallingInterval.id <= i.id + n_right)

    return [(x.chrom, x.start, x.end) for x in itertools.chain(left_intervals, [i], right_intervals)]


# define database model for storing the calling intervals
class ExacCallingInterval(peewee.Model):
    chrom = peewee.CharField(max_length=5, null=False, index=True)  # without 'chr'
    start = peewee.IntegerField()
    end = peewee.IntegerField()
    strand = peewee.CharField(max_length=1)
    name = peewee.CharField(max_length=500)

    class Meta:
        database = calling_intervals_db
        #primary_key = peewee.CompositeKey('chrom', 'start', 'end', 'strand', 'name')

indexes = (
    (('chrom', 'start', 'end', 'strand'), True),  # True means unique index
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

logging.info("db: %s(%s), autocommit=%s" % (type(calling_intervals_db).__name__, calling_intervals_db.database, calling_intervals_db.get_autocommit()))
logging.info("indexes: %s" % (calling_intervals_db.get_indexes("exaccallinginterval"),))

parse_exac_calling_intervals(args.exac_calling_intervals)
