"""
Utility module for parsing the exac info table into memory.

After this module is imported, the following global vars can be used to access
the data:

 EXAC_SAMPLE_ID_TO_BAM_PATH
 EXAC_SAMPLE_ID_TO_GVCF_PATH
 EXAC_SAMPLE_ID_TO_INCLUDE_STATUS

"""

import os

def parse_exac_info_table(info_table_path):
    # parse the ExAC info table to populate the following 3 dictionaries
    sample_id_to_bam_path = {}     # vcf sample_id => bam_path
    sample_id_to_gvcf_path = {}    # vcf sample_id => original GVCF produced as part of the ExAC pipeline run
    sample_id_include_status = {}  # vcf sample_id => True / False - based on "include" status in last column of info table

    with open(info_table_path) as info_table_file:
        header = next(info_table_file)
        for line in info_table_file:
            # info table has columns: vcf_sampleID, sampleID, ProjectID, ProjectName, Consortium, gvcf, bam, Include
            fields = line.strip('\n').split('\t')
            vcf_sample_id = fields[0]
            include_status = (fields[-1] == "YES")  # convert "YES" to boolean
            bam_path = fields[-2]
            gvcf_path = fields[-3]

            assert vcf_sample_id not in sample_id_to_bam_path, "duplicate sample id: %s" % vcf_sample_id

            sample_id_to_bam_path[vcf_sample_id] = bam_path
            sample_id_to_gvcf_path[vcf_sample_id] = gvcf_path
            sample_id_include_status[vcf_sample_id] = include_status

    return sample_id_to_bam_path, sample_id_to_gvcf_path, sample_id_include_status


import configargparse

p = configargparse.getArgumentParser()
p.add("--exac-info-table", help="Path of ExAC info table",
      default='/humgen/atgu1/fs03/lek/resources/ExAC/ExAC.r0.3_meta_Final.tsv')
args = p.parse_args()

assert os.path.isfile(args.exac_info_table), \
    "Couldn't find exac info table: %s" % args.exac_info_table

(EXAC_SAMPLE_ID_TO_BAM_PATH,
 EXAC_SAMPLE_ID_TO_GVCF_PATH,
 EXAC_SAMPLE_ID_TO_INCLUDE_STATUS) = parse_exac_info_table(args.exac_info_table)


assert len(EXAC_SAMPLE_ID_TO_BAM_PATH) == len(EXAC_SAMPLE_ID_TO_GVCF_PATH)
assert len(EXAC_SAMPLE_ID_TO_GVCF_PATH) == len(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS)
