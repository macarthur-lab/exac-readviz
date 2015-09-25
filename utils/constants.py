
# how many samples to show per het or hom-alt variant in the exac browser.
MAX_SAMPLES_TO_SHOW_PER_VARIANT = 5

# in ExAC v3 the max ref allele is 325bp long, and max alt allele is 365bp long, so use 500 to be safe.
MAX_ALLELE_SIZE = 500  # nucleotides

# in the ExAC v3 info table, the longest sample id is 36 chars long
MAX_VCF_SAMPLE_ID_SIZE = 50

# how many subdirectories to use to store the reassembled bams
NUM_OUTPUT_DIRECTORIES_L1 = 100
NUM_OUTPUT_DIRECTORIES_L2 = 10000

INCLUDE_N_ADJACENT_CALLING_REGIONS = 2

EXAC_CALLING_INTERVALS_PATH = "/seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list"
EXAC_INFO_TABLE_PATH = "/humgen/atgu1/fs03/lek/resources/ExAC/ExAC.r0.3_meta_Final.tsv"
EXAC_FULL_VCF_PATH = "/humgen/atgu1/fs03/konradk/exac/gqt/exac_all.vcf.gz"

GENCODE_BED_PATH = "/home/unix/weisburd/code/exac_readviz_scripts/data/gencode.v19.sorted.bed"
EXAC_CALLING_INTERVALS_BED_PATH = "/home/unix/weisburd/code/exac_readviz_scripts/data/exome_calling_regions.v1.bed"
SELF_CHAIN_BED_PATH = "/home/unix/weisburd/code/exac_readviz_scripts/data/self_chain.sorted.bed"

IGV_JAR_PATH = "/home/unix/weisburd/bin/IGV_2.3.44/igv.jar"
IGV_SCREEN_WIDTH = 300
IGV_TRACK_HEIGHT = 500
