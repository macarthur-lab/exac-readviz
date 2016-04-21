import pysam
from collections import defaultdict
from glob import glob
import re
import sys

total = affected = 0
for path in glob("../test/data/19/combined_*.bam"):
    f = pysam.AlignmentFile(path)

    pos_id_to_all_reads_counter = defaultdict(int)
    pos_id_to_interesting_reads_counter = defaultdict(int)
    #rg_to_reads = defaultdict(list)
    for r in f:
        if "D" not in r.cigarstring:
            continue

        tags = dict(r.tags)
        rg = tags['RG']  # eg. 19-51729103-CCCGG-C_hom2

        chrom, pos, ref, alt_and_g = rg.split("-")
        if len(ref) < 2:
            continue

        alt, hom_or_het_or_hemi = alt_and_g.split("_")
        pos_id = "%s-%s-%s-%s" % (chrom, pos, ref, alt)
        pos_id_to_all_reads_counter[pos_id] += 1

        pos = int(pos)
        if r.reference_start > pos or r.reference_end < pos:
            continue

        if not re.search("[^0-9]%dD" % (len(ref)-1), r.cigarstring):
            continue

        overlap = r.get_overlap(pos+1, pos+len(ref)-1)
        if not overlap: 
            continue   # deletion is exactly where it's expected

        if overlap == len(ref) - 1:
            continue   # deletion is in some other part of the read

        #actual_positions = set(r.get_reference_positions())
        #expected_deleted_positions = set(range(pos+1, pos+len(ref)))
        #print(rg + " deletion " + r.cigarstring + " " + str(r.get_overlap(pos+1, pos+len(ref)-1)) + " " + str(expected_deleted_positions) + " overlaps: " + str(actual_positions & expected_deleted_positions))

        pos_id_to_interesting_reads_counter[pos_id] += 1

        #rg_to_reads[rg].append(r)

    for pos_id in pos_id_to_all_reads_counter:
        total += 1
        if pos_id_to_interesting_reads_counter[pos_id] < 0.1*pos_id_to_all_reads_counter[pos_id]:
            continue

        print("http://exac.broadinstitute.org/variant/%s" % pos_id)
        affected += 1

    sys.stderr.write("%s out of %s (%0.1f)  %s\n" % (affected, total, 100.0*affected/total, pos_id))

