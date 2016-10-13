import sys
from tqdm import tqdm
import gzip
from utils.database import Variant, Sample
from utils.constants import MAX_ALLELE_SIZE

"""
Example:
1       12198   G       C       het     PROMIS20106530  /seq/picard_aggregation/C1996/PROMIS20106530/v1/PROMIS20106530.bam      /seq/sample_vcf/C1996/Exome/Homo_sapiens_assembly19/PROMIS20106530.20160617_0939.vcf.gz 0       27      10
"""

BUFFER_SIZE=1250

is_sample_id_on_fargo = {}
kristens_file = open("/humgen/atgu1/fs03/kristen/exac/metadata/status_fargo_all.txt")
header = next(kristens_file)  
print(header)
for line in kristens_file:
    fields = line.strip('\n').split('\t')
    sample_id = fields[1]
    is_sample_id_on_fargo[sample_id.lower()] = 1 if fields[-1].lower() == 'fargo' else 0
    #print(sample_id, fields[-1].lower(), is_sample_id_on_fargo[sample_id])

header = ['chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi', 'sample_id', 'original_bam_path', 'original_gvcf_path', 'sample_i', 'gq', 'dp']

filename = sys.argv[1]
print("Loading %s" % filename)

seen_variants = set()
variant_insert_buffer = []
sample_insert_buffer = []
for line in tqdm(gzip.open(filename), unit='row'):
    line = line.strip('\n')
    if not line:
        continue

    fields = line.split('\t')
    previous_variant = None
    
    try:
        assert len(header) == len(fields), "len(header) != len(fields)"

        row = dict(zip(header, fields))
        row['ref'] = row['ref'][:MAX_ALLELE_SIZE]
        row['alt'] = row['alt'][:MAX_ALLELE_SIZE]
        row['variant_id'] = "%s-%s-%s-%s" % ( row['chrom'], row['pos'], row['ref'], row['alt']) 
        row['priority'] = 0 if is_sample_id_on_fargo.get(row['sample_id'].lower()) else 1
        
        variant_uniq_id = row['variant_id'] + row['het_or_hom_or_hemi']
        if variant_uniq_id not in seen_variants:
            variant_insert_buffer.append({k: row[k] for k in ['chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi', 'variant_id']})
            seen_variants.add(variant_uniq_id)

        # sample
        assert int(row['dp']) >= 10, "int(row['dp']) < 10"
        assert int(row['gq']) >= 20, "int(row['gq']) < 20"
        
        del row['dp']
        del row['gq']

        sample_insert_buffer.append(row)

        if len(sample_insert_buffer) > BUFFER_SIZE:
            #with Sample._meta.database.atomic():
            Sample.insert_many(sample_insert_buffer).execute()
            sample_insert_buffer = []

        if len(variant_insert_buffer) > BUFFER_SIZE:
            #with Variant._meta.database.atomic():
            Variant.insert_many(variant_insert_buffer).execute()
            variant_insert_buffer = []


    except Exception as e:
        print("ERROR: %s on line %s" % (e, fields)) 
        raise

Variant.insert_many(variant_insert_buffer).execute()
Sample.insert_many(sample_insert_buffer).execute()
    

