import sys
from tqdm import tqdm
import gzip
from collections import defaultdict
from utils.database import Sample
from utils.constants import MAX_ALLELE_SIZE
from utils.minimal_representation import get_minimal_representation
from utils.exac_info_table import compute_latest_bam_path
from utils.file_utils import does_file_exist

"""
Example:
1       12198   G       C       het     PROMIS20106530  /seq/picard_aggregation/C1996/PROMIS20106530/v1/PROMIS20106530.bam      /seq/sample_vcf/C1996/Exome/Homo_sapiens_assembly19/PROMIS20106530.20160617_0939.vcf.gz 0       27      10
"""

BUFFER_SIZE=1250

class SampleGenomes(Sample):
    class Meta:
        db_table = 'sample_genomes'

class SampleExomes(Sample):
    class Meta:
        db_table = 'sample_exomes'

#filename = "/humgen/atgu1/fs03/weisburd/exac_readviz_data_v2/laurents_tables_v2/gnomad.readviz.txt.bgz"
#exome_or_genome = "G"
#Sample = SampleGenomes

filename = "/humgen/atgu1/fs03/weisburd/exac_readviz_data_v2/laurents_tables_v2/exacv2.readviz.txt.bgz"
exome_or_genome = "E"
Sample = SampleExomes

print("Loading %s" % filename)
print("Table name: " + Sample._meta.db_table)
#sys.exit(0)

header = ['chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi', 'sample_id', 'original_bam_path', 'original_gvcf_path', 
          'project_id', 'project_description', 'sex', 'population', 'sample_i', 'gq', 'dp', 'dp_ref', 'dp_alt']


file_exists = {}
sample_insert_buffer = []
counter = defaultdict(int)
f = gzip.open(filename)
print("Skipping header line: " + next(f))
for line in tqdm(f, unit='row'):
    line = line.strip('\n')
    if not line:
        continue

    counter['total'] += 1

    fields = line.split('\t')

    try:
        assert len(header) == len(fields), "len(header) != len(fields)"

        row = dict(zip(header, fields))
        row['exome_or_genome'] = exome_or_genome

        int_pos = int(row['pos'])
        int_pos, row['ref'], row['alt'] = get_minimal_representation(int_pos, row['ref'], row['alt'])
        row['ref'] = row['ref'][:MAX_ALLELE_SIZE]
        row['alt'] = row['alt'][:MAX_ALLELE_SIZE]
        row['variant_id'] = "%s-%s-%s-%s" % ( row['chrom'], row['pos'], row['ref'], row['alt']) 

        if len(row['ref']) > 1 and len(row['alt']) > 1:
            print("ERROR: couldn't min-rep " + row['variant_id'])


        row['priority'] = 0 
        row['original_bam_path'] = compute_latest_bam_path(row['original_bam_path'])

        if row['original_bam_path'] not in file_exists:
            file_exists[ row['original_bam_path'] ] = int(does_file_exist(row['original_bam_path']))
            counter['files_that_exist'] = sum(file_exists.values())
            counter['total_files'] = len(file_exists)
            counter['files_that_are_missing'] = counter['total_files'] - counter['files_that_exist']
            if not  file_exists[ row['original_bam_path'] ]:
                print("\nMissing bam: " + row['original_bam_path'])
            
        row['pos'] = str(int_pos)
        row['pos_mod_1000'] = int_pos % 1000        
        row['sex'] = row['sex'][0].upper()

        # sample
        assert int(row['dp']) >= 10, "int(row['dp']) < 10"
        assert int(row['gq']) >= 20, "int(row['gq']) < 20"
        
        sample_insert_buffer.append(row)
        if len(sample_insert_buffer) > BUFFER_SIZE:
            sql = Sample.insert_many(sample_insert_buffer)
            sql.execute()
            sample_insert_buffer = []
            print("\n" + str(dict(counter)))

    except Exception as e:
        print("ERROR: %s on line %s" % (e, fields)) 
        raise

Sample.insert_many(sample_insert_buffer).execute()
