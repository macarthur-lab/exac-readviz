"""
Utility module for parsing the exac info table into memory.

After this module is imported, the following global vars can be used to access
the data:
 TCGA_SAMPLE_ID_TO_BAM_PATH
 TCGA_SAMPLE_ID_TO_GVCF_PATH

 EXAC_SAMPLE_ID_TO_BAM_PATH
 EXAC_SAMPLE_ID_TO_GVCF_PATH
 EXAC_SAMPLE_ID_TO_INCLUDE_STATUS
 
"""

import logging
import os
import re
from tqdm import tqdm

from utils.constants import TCGA_NEW_BAM_PATHS

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


def parse_exac_pop_sex_table(pop_sex_table_path):
    sample_id_to_population = {}     # vcf sample_id => population (eg. 'AFR')
    sample_id_to_sex = {}     # vcf sample_id => 'm' or 'f'
    with open(pop_sex_table_path) as pop_sex_table_file:
        header = ['sample_id', 'pop', 'sex']
        for i, line in enumerate(pop_sex_table_file):
            fields = line.strip('\n').split('\t')
            vcf_sample_id = fields[0]

            pop = fields[1]
            sample_id_to_population[vcf_sample_id] = pop

            sex = fields[2]
            if sex == 'Male':
                sex = 'm'
            elif sex == 'Female':
                sex = 'f'
            else:
                raise ValueError("line %s has unexpected value for sex: '%s'" % (i, sex))

            sample_id_to_sex[vcf_sample_id] = sex

    return sample_id_to_population, sample_id_to_sex
    

try:
    from utils.constants import EXAC_INFO_TABLE_PATH, EXAC_POP_SEX_TABLE_PATH
    assert os.path.isfile(EXAC_INFO_TABLE_PATH), \
        "Couldn't find exac info table: %s" % EXAC_INFO_TABLE_PATH

    assert os.path.isfile(EXAC_POP_SEX_TABLE_PATH), \
        "Couldn't find exac pop sex table: %s" % EXAC_POP_SEX_TABLE_PATH
    

    (EXAC_SAMPLE_ID_TO_BAM_PATH,
     EXAC_SAMPLE_ID_TO_GVCF_PATH,
     EXAC_SAMPLE_ID_TO_INCLUDE_STATUS) = parse_exac_info_table(EXAC_INFO_TABLE_PATH)

    assert len(EXAC_SAMPLE_ID_TO_BAM_PATH) == len(EXAC_SAMPLE_ID_TO_GVCF_PATH)
    assert len(EXAC_SAMPLE_ID_TO_GVCF_PATH) == len(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS)
    

    n_total = len(EXAC_SAMPLE_ID_TO_BAM_PATH)
    n_include_true = sum(EXAC_SAMPLE_ID_TO_INCLUDE_STATUS.values())
    
    logging.info("Loaded %s" % EXAC_INFO_TABLE_PATH)
    logging.info("INCLUDE_STATUS = True in %d out of %d (%0.1f%%) samples" % (
        n_include_true, n_total, 100*n_include_true/float(n_total)))


    (EXAC_SAMPLE_ID_TO_POP,
     EXAC_SAMPLE_ID_TO_SEX) = parse_exac_pop_sex_table(EXAC_POP_SEX_TABLE_PATH)

    n_male = len([s for s in EXAC_SAMPLE_ID_TO_SEX.values() if s == 'm'])
    n_female = len([s for s in EXAC_SAMPLE_ID_TO_SEX.values() if s == 'f'])

    assert n_male + n_female == len(EXAC_SAMPLE_ID_TO_POP), \
        "n_male (%s) + n_female (%s) != len(EXAC_SAMPLE_ID_TO_POP) (%s)" % (n_male, n_female, len(EXAC_SAMPLE_ID_TO_POP))

    logging.info("Loaded %s" % EXAC_POP_SEX_TABLE_PATH)
    logging.info("%d male, %d female" % (n_male, n_female))
except Exception, e:
    print("WARNING: " + str(e))

TCGA_SAMPLE_ID_TO_BAM_PATH = TCGA_SAMPLE_ID_TO_GVCF_PATH = {}
try:
    # Parse the new TCGA bam paths
    with open(TCGA_NEW_BAM_PATHS) as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            tcga_sample_id = fields[0]

            tcga_gvcf_path = fields[5]
            tcga_bam_path = fields[-1]


            #if not os.path.isfile(tcga_gvcf_path): print("ERROR: vcf file not found: " + tcga_gvcf_path)
            #if not os.path.isfile(tcga_bam_path): print("ERROR: bam file not found: " + tcga_gvcf_path)

            TCGA_SAMPLE_ID_TO_BAM_PATH[tcga_sample_id] = tcga_bam_path
            TCGA_SAMPLE_ID_TO_GVCF_PATH[tcga_sample_id] = tcga_gvcf_path
except Exception, e:
    print("WARNING: " + str(e))

    
try:
    for tcga_sample_id in TCGA_SAMPLE_ID_TO_BAM_PATH:
        assert tcga_sample_id in EXAC_SAMPLE_ID_TO_BAM_PATH
        assert tcga_sample_id in EXAC_SAMPLE_ID_TO_GVCF_PATH

        EXAC_SAMPLE_ID_TO_BAM_PATH[tcga_sample_id] = TCGA_SAMPLE_ID_TO_BAM_PATH[tcga_sample_id]
        assert EXAC_SAMPLE_ID_TO_GVCF_PATH[tcga_sample_id] == tcga_gvcf_path

        print('---\n%s\n%s' % (EXAC_SAMPLE_ID_TO_BAM_PATH[tcga_sample_id], tcga_bam_path))

except Exception, e:
    print("WARNING: " + str(e))
            

BAM_PATH_REGEXP = re.compile("/v[0-9]{1,2}/")

def compute_latest_bam_path(bam_path):
    """Compute current bam path"""

    #if "/cga/pancan2/picard_bams/ext_tcga" in bam_path:
    #    bam_path = bam_path.replace("/cga/pancan2/", "/cga/fh/cga_pancan2/")

    if "/C1437/" in bam_path:
        d, n = os.path.split(bam_path)
        prefix = n.replace(".bam", "")
        new_prefix = prefix.replace("_", "")
        bam_path = os.path.join(d.replace(prefix, new_prefix), new_prefix+".bam")
    elif "CONT_" in bam_path:
        bam_path = bam_path.replace("CONT_", "CONT")

    bam_path = BAM_PATH_REGEXP.sub("/current/", bam_path)  # get latest version of the bam

    return bam_path


def lookup_original_bam_path(sample_id):
    """Look up the bam path for the given sample id.

    Args:
      sample_id: vcf sample id
    Return:
      The original bam path
    """

    # work-arounds for relocated .bams taken from @birndle's igv_spot_checking script
    bam_path = EXAC_SAMPLE_ID_TO_BAM_PATH[sample_id]

    return compute_latest_bam_path(bam_path)
