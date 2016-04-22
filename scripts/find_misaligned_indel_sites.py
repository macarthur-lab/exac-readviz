import pysam
from collections import defaultdict
from glob import glob
import os
import re
import sys
import peewee

from utils.database import _SharedMeta
from utils.constants import MAX_ALLELE_SIZE


# get the sample
class _SharedVariantFields(_SharedMeta):
    chrom = peewee.CharField(max_length=5, null=False, index=True)
    pos = peewee.IntegerField()
    ref= peewee.CharField(max_length=MAX_ALLELE_SIZE)
    alt = peewee.CharField(max_length=MAX_ALLELE_SIZE)
    het_or_hom = peewee.CharField(max_length=4)
    #het_or_hom_or_hemi = peewee.CharField(max_length=4)


class Sample(_SharedVariantFields):
    sample_id = peewee.CharField()
    sample_i = peewee.IntegerField(null=True)
    hc_succeeded = peewee.BooleanField(default=0, index=True)
    original_bam_path = peewee.TextField(null=True)

    def toJSON(self):
        return {'chrom' : self.chrom, 'pos': self.pos, 'ref': self.ref, 'alt': self.alt, 'het_or_hom' : self.het_or_hom, 'sample_i' : self.sample_i}

    def __hash__(self):
        return hash(str(self.toJSON()))


    class Meta:
        db_table = 'previous_2016_03_27__sample'


print("\t".join(["pos_id", "sample_id", "bam_path"]))

total = affected = 0
for path in glob(os.path.join(os.path.dirname(__file__), "../test/data/19/combined_*.bam")):
    f = pysam.AlignmentFile(path)

    pos_id_to_all_reads_counter = defaultdict(int)
    pos_id_to_interesting_reads_counter = defaultdict(int)
    pos_id_to_sample_records = defaultdict(set)
    for r in f:
        if "D" not in r.cigarstring:
            continue

        tags = dict(r.tags)
        rg = tags['RG']  # eg. 19-51729103-CCCGG-C_hom2

        chrom, pos, ref, alt_and_g = rg.split("-")
        if len(ref) < 2:
            continue

        alt, het_or_hom_or_hemi_and_i = alt_and_g.split("_")
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


        # get the original bam paths
        het_or_hom_or_hemi = het_or_hom_or_hemi_and_i[0:-1]
        sample_i = int(het_or_hom_or_hemi_and_i[-1])

        try:
            s = Sample.get(chrom=chrom, pos=pos, ref=ref, alt=alt, het_or_hom=het_or_hom_or_hemi, sample_i=sample_i, hc_succeeded=1)
            pos_id_to_sample_records[pos_id].add(s)
        except Exception, e:
            sys.stderr.write("Exception: %s \n" % e)
            continue

    for pos_id in pos_id_to_all_reads_counter:
        total += 1
        if pos_id_to_interesting_reads_counter[pos_id] < 0.1*pos_id_to_all_reads_counter[pos_id]:
            continue

        #print(rg)
        for s in pos_id_to_sample_records[pos_id]:
            url = "http://exac.broadinstitute.org/variant/%s" % pos_id
            print("\t".join([pos_id, s.sample_id, s.original_bam_path]))

        affected += 1

    sys.stderr.write("%s out of %s (%0.1f)  %s\n" % (affected, total, 100.0*affected/total, pos_id))

