import logging
import os

from utils.file_utils import does_file_exist

DB_HOST = 'exac-dev'
DB_PORT = 3307 
DB_USER = 'root'

DATA_DIR_PREFIX = '/humgen/atgu1/fs03/weisburd/exac_readviz_scripts_data' 
BIN_DIR_PREFIX = os.path.join(DATA_DIR_PREFIX, 'bin') 
BAM_OUTPUT_DIR = "/broad/hptmp/exac_readviz_backend/"


# how many samples to show per het, hom-alt or hemizygous variant in the exac browser.
MAX_SAMPLES_TO_SHOW_PER_VARIANT = 5

# in ExAC v3 the max ref allele is 325bp long, and max alt allele is 365bp long, so use 500 to be safe.
MAX_ALLELE_SIZE = 500  # nucleotides

# in the ExAC v3 info table, the longest sample id is 36 chars long
MAX_VCF_SAMPLE_ID_SIZE = 50

# how many subdirectories to use to store the reassembled bams
#NUM_OUTPUT_DIRECTORIES_L1 = 100
#NUM_OUTPUT_DIRECTORIES_L2 = 10000
NUM_OUTPUT_DIRECTORIES_L1 = 1000

INCLUDE_N_ADJACENT_CALLING_REGIONS = 2

EXAC_CALLING_INTERVALS_PATH = os.path.join(DATA_DIR_PREFIX, "exome_calling_regions.v1.interval_list")
EXAC_INFO_TABLE_PATH = os.path.join(DATA_DIR_PREFIX, "ExAC.r0.3_meta_Final.tsv")
EXAC_POP_SEX_TABLE_PATH = os.path.join(DATA_DIR_PREFIX, "samples_pop_sex.tsv")
EXAC_FULL_VCF_PATH = os.path.join(DATA_DIR_PREFIX, "exac_all.vcf.gz")
EXAC_SITES_VCF_PATH = os.path.join(DATA_DIR_PREFIX, "ExAC.r0.3.1.sites.vep.vcf.gz")
GENCODE_EXAC_GTF_PATH = os.path.join(DATA_DIR_PREFIX, "gencode.gtf.gz")

PICARD_JAR_PATH = os.path.join(BIN_DIR_PREFIX,"picard.jar")  # used for sorting bam

GATK_JAR_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)),
#    "bin/noMQ0sInBamout/GenomeAnalysisTK.jar"))
    "bin/GATK_noMQ0sInBamout_fixRealign.jar"))

TCGA_NEW_BAM_PATHS = os.path.join(DATA_DIR_PREFIX, "TCGA_external_cghublink_batch1.tsv")

# used for igv screenshots
GENCODE_BED_PATH = os.path.join(DATA_DIR_PREFIX, "gencode.v19.sorted.bed")
EXAC_CALLING_INTERVALS_BED_PATH = os.path.join(DATA_DIR_PREFIX, "exome_calling_regions.v1.bed")
SELF_CHAIN_BED_PATH = os.path.join(DATA_DIR_PREFIX, "self_chain.sorted.bed")

IGV_JAR_PATH = os.path.join(BIN_DIR_PREFIX, "IGV_2.3.44/igv.jar")
IGV_SCREEN_WIDTH = 300
IGV_TRACK_HEIGHT = 500

all_files_exist = True
for path in (EXAC_CALLING_INTERVALS_PATH, EXAC_INFO_TABLE_PATH,
             EXAC_POP_SEX_TABLE_PATH, EXAC_FULL_VCF_PATH, EXAC_SITES_VCF_PATH,
             EXAC_SITES_VCF_PATH, GENCODE_EXAC_GTF_PATH,
             PICARD_JAR_PATH, GATK_JAR_PATH,
             TCGA_NEW_BAM_PATHS):
    if not does_file_exist(path, use_cache=False):
        logging.error("ERROR: file not found: " + path)
        all_files_exist = False

if not all_files_exist:
    raise Exception("Critical files not found")
