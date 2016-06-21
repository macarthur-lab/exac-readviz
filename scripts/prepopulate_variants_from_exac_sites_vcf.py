import argparse
import collections
import gzip
import logging
from peewee import fn
import pysam

from utils.constants import EXAC_SITES_VCF_PATH, MAX_SAMPLES_TO_SHOW_PER_VARIANT
from utils.database import Variant
from utils.exac_vcf import create_vcf_row_parser, create_variant_iterator_from_vcf

p = argparse.ArgumentParser()
p.add_argument("-p", "--exac-sites-vcf-path", help="ExAC sites VCF path", default=EXAC_SITES_VCF_PATH)
p.add_argument("--chrom", help="Chromosome to process")
p.add_argument("--restart", action="store_true", help="Restart loading from where it left off")

args = p.parse_args()

Variant.create_table(fail_silently=True)

tabix_file = pysam.TabixFile(filename=args.exac_sites_vcf_path, parser=pysam.asTuple())
last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
if args.chrom:
    start_at_pos = 1
    if args.restart:
        variants_with_max_pos = list(Variant.select(fn.Max(Variant.pos).alias('pos')).where(Variant.chrom==args.chrom).limit(1))
        if variants_with_max_pos:
            start_at_pos = variants_with_max_pos[0].pos

    logging.info("Processing chromosome: %s starting at %s" % (args.chrom, start_at_pos))
    vcf_row_iterator = tabix_file.fetch(args.chrom, start_at_pos, 5*10**8)
else:
    vcf_row_iterator = pysam.tabix_iterator(gzip.open(args.exac_sites_vcf_path), parser=pysam.asTuple())
    logging.info("Processing all chromosomes")


# create variant iterator
parse_vcf_row = create_vcf_row_parser(last_header_line)

variant_iterator = create_variant_iterator_from_vcf(vcf_row_iterator, parse_vcf_row)

# iterate over variants
counters = collections.defaultdict(int)
for chrom, pos, ref, alt, alt_allele_index, het_or_hom_or_hemi, n_expected_samples, all_genotypes_in_row in variant_iterator:

    counters["all_alleles"] += 1

    # check if allele has been processed already (skip if yes)
    vr, created = Variant.get_or_create(chrom=chrom, pos=pos, ref=ref, alt=alt, het_or_hom_or_hemi=het_or_hom_or_hemi)

    if created:
        vr.variant_id = "%s-%s-%s-%s" % (chrom, pos, ref, alt)
        vr.finished = 0
        vr.n_available_samples = 0
        vr.n_expected_samples = min(n_expected_samples, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
        vr.save()

        counters["created_alleles"] += 1
        if counters["created_alleles"] % 100 == 0:
            stats = ", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])])
            logging.info("created %s-%s-%s-%s %s  n_expected=%s  %s" % (vr.chrom, vr.pos, vr.ref, vr.alt, vr.het_or_hom_or_hemi, vr.n_expected_samples, stats))

    else:
        if vr.finished and vr.n_expected_samples != n_expected_samples:
            logging.error("%s: variant record n_expected_samples == %d, while n_expected_samples in ExAC sites vcf == %d, %s" % (
                vr.het_or_hom_or_hemi, vr.n_expected_samples, n_expected_samples, all_genotypes_in_row))

print("Finished")
