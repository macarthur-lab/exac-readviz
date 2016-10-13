"""Remove samples from sample table that were dropped due to QC, etc."""

import os
from tqdm import tqdm
from utils.exac_info_table import compute_latest_bam_path
from utils.database import Sample, _readviz_db


#all_current_samples = set(s[0] for s in _readviz_db.execute_sql("select distinct sample_id from sample"))
#print(all_current_samples)
#print(len(all_current_samples))

#sys.exit(0)
test_run = False
for line in open(os.path.join(os.path.dirname(__file__), "samples_to_drop.txt")):
    sample_id = line.strip()        
    print("Processing %s" % sample_id)
    samples_to_delete = list(Sample.select().where(Sample.sample_id == sample_id).order_by(Sample.sample_id.asc()).execute())
    for s in samples_to_delete:
        samples_to_adjust = list(Sample.select().where(
                Sample.chrom == s.chrom, 
                Sample.pos == s.pos, 
                Sample.ref == s.ref, 
                Sample.alt == s.alt, 
                Sample.het_or_hom_or_hemi == s.het_or_hom_or_hemi, 
                Sample.sample_i > s.sample_i).order_by(Sample.sample_i.asc()))

        if len(list(Sample.select().where(
                Sample.chrom == s.chrom,
                Sample.pos == s.pos,
                Sample.ref == s.ref,
                Sample.alt == s.alt,
                Sample.het_or_hom_or_hemi == s.het_or_hom_or_hemi,
                Sample.sample_i == s.sample_i))) > 1:
            if not test_run:
                Sample.delete().where(Sample.chrom == s.chrom,
                                      Sample.pos == s.pos,
                                      Sample.ref == s.ref,
                                      Sample.alt == s.alt,
                                      Sample.het_or_hom_or_hemi == s.het_or_hom_or_hemi, 
                                      Sample.sample_i == s.sample_i, 
                                      Sample.sample_id == s.sample_id).execute()
                print("    deleted %s %s sample%s, skipping %s previously-updated samples" % (s.variant_id, s.het_or_hom_or_hemi, s.sample_i, len(samples_to_adjust)))
                continue

        for s2 in list(samples_to_adjust):
            s2.priority = 0
            s2.started = 0
            s2.started_time = None
            s2.finished = 0
            s2.finished_time = None
            s2.comments = None
            s2.username = None
            s2.sample_i = s2.sample_i - 1
            s2.hc_error_code = None
            s2.hc_error_text = None
            s2.hc_command_line = None
            s2.hc_succeeded = 0
            if not test_run:
                s2.save()
        print("    deleted %s %s sample%s, updated %s samples" % (s.variant_id, s.het_or_hom_or_hemi, s.sample_i, len(samples_to_adjust)))
        if not test_run:
            Sample.delete().where(Sample.chrom == s.chrom,
                Sample.pos == s.pos,
                Sample.ref == s.ref,
                Sample.alt == s.alt,
                Sample.het_or_hom_or_hemi == s.het_or_hom_or_hemi, 
                Sample.sample_i == s.sample_i, 
                Sample.sample_id == s.sample_id).execute()


