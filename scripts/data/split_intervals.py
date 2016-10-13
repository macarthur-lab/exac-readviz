#!/usr/bin/env python2.7

import argparse
import datetime
import os
import peewee
import getpass
import random
import signal
import slugify
import subprocess
from collections import defaultdict

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

p = argparse.ArgumentParser()
p.add_argument("-isize", "--interval-size", help="max interval size", type=int, required=True)
p.add_argument("interval_list", help="input intervals file to split")
args, unknown_args = p.parse_known_args()

logging.info("args: %s" % str(args))

output_filename = str(args.interval_size) + "bp_" + os.path.basename(args.interval_list)
print(output_filename)

# create intervals
intervals = []

input_bases = 0
output_bases = 0
with open(str(args.interval_list)) as interval_list_file, open(output_filename, "w") as f_out:
    for line in interval_list_file:
        if line.startswith("@"):
            f_out.write(line)
            continue

        fields = line.strip('\n').split("\t")
        chrom, start, end, strand, name = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
        input_bases += end - start + 1
        for prev_interval in intervals[-4:]:
            prev_start = prev_interval[1]
            prev_end = prev_interval[2]
            if end >= prev_start and start <= prev_end:
                print("%s:%s-%s overlaps prev interval %s:%s-%s" % (chrom, start, end, chrom, prev_start, prev_end))
        intervals.append((chrom, start, end, strand, name))
        
    logging.info("Parsed %s intervals from %s" % (len(intervals), args.interval_list))

    # split intervals so they are no bigger than args.interval_size
    final_intervals = []
    for chrom, start, end, strand, name in intervals:
        if end - start > args.interval_size:
            #print("Breaking up large interval %s:%s-%s" % (chrom, start, end))
            while end - start > args.interval_size:
                #print("   %s:%s-%s" % (chrom, start, start + args.interval_size-1))
                new_end = start + args.interval_size-1
                f_out.write("\t".join([chrom, str(start), str(new_end), strand, name]) + "\n")
                output_bases += new_end - start + 1

                final_intervals.append({"chrom": chrom, "start_pos": start, "end_pos": new_end})
                start = start + args.interval_size

        
        final_intervals.append({"chrom": chrom, "start_pos": start, "end_pos": end})
        f_out.write("\t".join([chrom, str(start), str(end), strand, name]) + "\n")
        output_bases += end - start + 1

    logging.info("Broke the %s intervals into %s intervals of size %s or less" % (len(intervals), len(final_intervals), args.interval_size))

print("input  bases: %s" % input_bases)
print("output bases: %s" % output_bases)
