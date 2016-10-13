import os
from tqdm import tqdm
from utils.exac_info_table import compute_latest_bam_path
from utils.database import Sample

input_file = "all_bams_paths.txt"
output_file = "missing_bams.txt"

with open(output_file, "w") as f:
    counter = 0
    if os.path.isfile(input_file):
        all_bam_paths = [line.strip('\n') for line in open(input_file)]
    else:
        all_bam_paths = [s.original_bam_path for s in tqdm(Sample.select(Sample.original_bam_path).group_by(Sample.original_bam_path).execute())]
        with open(input_file, "w") as f:
            for path in all_bam_paths:
                f.write("%s\n" % path)

    with open(output_file, "w") as f:
        for bam_path in tqdm(all_bam_paths, unit=" bams"):
            bam_path = compute_latest_bam_path(bam_path)
            #if "external-data" in bam_path:
            #    continue

            if not os.path.isfile(bam_path):
                print("ERROR: bam not found: %s" % bam_path)
                f.write("%s\n" % bam_path)
                f.flush()
                counter += 1

    print("Done - %s out of %s bams missing" % (counter, len(all_bam_paths)))

