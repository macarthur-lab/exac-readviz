"""
This module parses the calling intervals file into a database table (if
this hasn't been done already). The table is indexed by (chrom, start, end).
"""
import configargparse
import logging
import os
import sys

from utils.constants import EXAC_CALLING_INTERVALS_PATH
from utils.database import ExacCallingInterval

# utility functions
def parse_exac_calling_intervals(exac_calling_intervals_path):
    """Parses the exac calling regions .intervals file into the
    'ExacCallingInterval' table.

    Args:
        exac_calling_intervals_path: A path like
        "/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list"
    """

    # check if already done
    with open(exac_calling_intervals_path) as f:
        lines_in_file = sum(1 for line in f if not line.startswith("@"))

    records_in_db = ExacCallingInterval.select().count()
    if records_in_db > 0:
        message = "%d out of %d intervals have been loaded from %s into %s('%s')" % (
            records_in_db, lines_in_file, exac_calling_intervals_path,
            type(ExacCallingInterval.database).__name__,
            ExacCallingInterval.database)
        if lines_in_file == records_in_db:
            logging.info("all " + message)
            _print_stats()
            return
        else:
            logging.info("Only " + message)
    else:
        logging.info("loading exac calling intervals from %s" % exac_calling_intervals_path)

    # parse exac_calling_intervals_file into the database table
    with open(exac_calling_intervals_path) as f:
        previous_chrom = None
        previous_start_pos = None
        line_counter = 0

        for line in f:
            # skip header
            if line.startswith("@"):
                continue

            line_counter += 1
            if line_counter % 10000 == 0:
                logging.info("%d lines processed" % line_counter)

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

    _print_stats()
    logging.info("done loading")

def _print_stats():
    for row in ExacCallingInterval.raw("select count(*) as n, min(end-start) as min_size, max(end-start) as max_size, avg(end-start) as avg_size from " + ExacCallingInterval.__name__.lower()):
        logging.info("stats: %s intervals (range: %s to %s, mean: %s)" % (row.n, row.min_size, row.max_size, int(row.avg_size)))


def _get_interval(chrom, pos):
    """Utility method - returns the ExacCallingInterval object that overlaps the
    given chrom, pos (1-based, inclusive)
    """
    intervals = list(
        ExacCallingInterval.select().where((
                ExacCallingInterval.chrom == chrom
            ) & (
                ExacCallingInterval.start <= pos
            ) & (
                ExacCallingInterval.end >= pos
            ))
    )

    if not intervals:
        raise ValueError(
            "Variant %s:%s is not overlapped by any calling intervals" % (chrom, pos))

    assert len(intervals) < 2, "Multiple calling intervals overlap variant %s:%s - %s" % (
        chrom, pos, intervals)

    return intervals[0]


def get_overlapping_calling_interval(chrom, pos):
    """Returns the for the calling interval that overlaps
    the given chrom, pos (1-based, inclusive).
    """
    return _get_interval(chrom, pos)



def get_adjacent_calling_intervals(chrom, pos, n_left=1, n_right=1):
    """Returns the calling interval that overlaps the given chrom: pos (1-based,
    inclusive), and also up to n_left intervals to the left of this interval
    (in order by chrom & start position), as well as up to n_right intervals
    to the right of the interval. If the overlapping interval is at the very
    beginning or very end of a chromosome, fewer than n_left or n_right
    intervals may be included.

    Args:
      chrom: variant chrom
      pos: variant pos
      n_left: returned list will include this many intervals to the left of the
        interval that overlaps the given variant
      n_right: returned list will include this many intervals to the right of the
        interval that overlaps the given variant

    Return:
        3-tuple (x,y,z) where
            x = up to n_left ExacCallingInterval objects representing neighboring intervals to the left
            y = the ExacCallingInterval that overlaps the given variant
            z = up to n_right ExacCallingInterval objects representing neighboring intervals to the right
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

    return list(left_intervals), i, list(right_intervals)


if __name__ == "__main__":

    # parse the exac calling intervals file into a database table
    assert os.path.isfile(EXAC_CALLING_INTERVALS_PATH), \
        "Couldn't find: %s" % EXAC_CALLING_INTERVALS_PATH

    parse_exac_calling_intervals(EXAC_CALLING_INTERVALS_PATH)
