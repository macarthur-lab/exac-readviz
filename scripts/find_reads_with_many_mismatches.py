import argparse
from collections import defaultdict
import os
import pysam
import sys


p = argparse.ArgumentParser()
p.add_argument("combined_bam", nargs="+")
args = p.parse_args()

print("%s files" % len(args.combined_bam))


fasta = pysam.Fastafile("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
total_variants = set()
variants = defaultdict(int)
mismatch_counts = []
for bam_path in args.combined_bam:
    f = pysam.AlignmentFile(bam_path)
    for r in f:
        tags = dict(r.tags)
        rg = tags['RG']  # eg. 19-51729103-CCCGG-C_hom2

        total_variants.add(rg)

        if len(r.cigarstring) > 4 or 'I' in r.cigarstring or 'D' in r.cigarstring or 'M' not in r.cigarstring:
            continue

        variant_pos = int(rg.split('-')[1])
        if not ((r.reference_start - 50 < variant_pos) and (r.reference_end + 50 > variant_pos)):
            continue


        read_sequence = r.query_alignment_sequence
        reference_sequence = fasta.fetch(r.reference_name, r.reference_start, r.reference_end)
        if len(read_sequence) != len(reference_sequence):
            print("ERROR: Different lengths - cigar: %s" % r.cigarstring)
            continue
        if read_sequence != reference_sequence:
            mismatches = [i for i, (a,b) in enumerate(zip(read_sequence, reference_sequence)) if a != b and a != 'N' and b != 'N']
            if len(mismatches) > 10:
                variants[rg] += 1
                mismatch_counts.append(1 - len(mismatches)/float(len(read_sequence)))
                print("#%s: %s %s - %s mismtaches - read: %s:%s %s" % (len(variants), rg.replace("_", " "), r.cigarstring, len(mismatches), r.reference_name, r.reference_start, r.query_name))

    class FakeOptions:
        def __getattr__(self, key):
            return None
    if len(set(variants.values())) > 1:
        print("Bad reads per variant:")
        os.system("echo '%s' | histogram.py" % '\n'.join(map(str, variants.values())))
    if len(set(mismatch_counts)) > 1:
        print("Fraction of read that matches the reference:")
        os.system("echo '%s' | histogram.py" % '\n'.join(map(str, mismatch_counts)))

    print("%s out of %s variants (%s)" % (len(variants), len(total_variants), len(variants)/float(len(total_variants))))

print("variants with mismatching reads :")
print("\n".join(variants))

