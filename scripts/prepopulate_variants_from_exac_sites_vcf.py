import argparse
import collections
import gzip
import logging
from peewee import fn
import pysam
import sys

from utils.constants import EXAC_SITES_VCF_PATH, MAX_SAMPLES_TO_SHOW_PER_VARIANT
from utils.database import Variant
from utils.exac_vcf import create_vcf_row_parser
from utils.minimal_representation import get_minimal_representation

p = argparse.ArgumentParser()
p.add_argument("-p", "--exac-sites-vcf-path", help="ExAC sites VCF path", default=EXAC_SITES_VCF_PATH)
p.add_argument("--chrom", help="Chromosome to process", required=True)
args = p.parse_args()

start_at_pos = 1
#for v in Variant.select(fn.Max(Variant.pos).alias('max_pos')).where(Variant.chrom==args.chrom):
#    start_at_pos = v.max_pos


tabix_file = pysam.TabixFile(filename=args.exac_sites_vcf_path, parser=pysam.asTuple())
last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
if args.chrom:
    logging.info("Processing chromosome: %s starting at %s" % (args.chrom, start_at_pos))
    vcf_iterator = tabix_file.fetch(args.chrom, start_at_pos, 5*10**8)
else:
    vcf_iterator = pysam.tabix_iterator(gzip.open(args.exac_sites_vcf_path), parser=pysam.asTuple())
    logging.info("Processing all chromosomes")


parse_vcf_row = create_vcf_row_parser(last_header_line)

counters = collections.defaultdict(int)
for row in vcf_iterator:
    chrom, pos, ref, alt_alleles, n_het_list, n_hom_list, n_hemi_list, all_genotypes_in_row = parse_vcf_row(row)
    counters["sites"] += 1

    # iterate over alt alleles (in case this row is multi-allelic)
    for alt_allele_index, (alt, n_het, n_hom, n_hemi) in enumerate(zip(alt_alleles, n_het_list, n_hom_list, n_hemi_list)):

        minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(pos, ref, alt)

        counters["all_alleles"] += 1

        # choose het and hom-alt samples to display for this allele
        for het_or_hom_or_hemi in ["het", "hom", "hemi"]:
            # check if allele has been processed already (skip if yes)
            vr, created = Variant.get_or_create(chrom=chrom, pos=minrep_pos, ref=minrep_ref, alt=minrep_alt, het_or_hom_or_hemi=het_or_hom_or_hemi)

            if created:
                vr.finished = 0
                vr.n_available_samples = 0
                if het_or_hom_or_hemi == "het":
                    vr.n_expected_samples = min(n_het, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
                elif het_or_hom_or_hemi == "hom":
                    vr.n_expected_samples = min(n_hom, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
                elif het_or_hom_or_hemi == "hemi":
                    vr.n_expected_samples = min(n_hemi, MAX_SAMPLES_TO_SHOW_PER_VARIANT)
                else:
                    raise ValueError("Unexpected value for het_or_hom_or_hemi: %s" % str(het_or_hom_or_hemi))

                vr.save()

                counters["created_alleles"] += 1
                if counters["created_alleles"] % 100 == 0:
                    logging.info("created %s-%s-%s-%s %s  %s" % (vr.chrom, vr.pos, vr.ref, vr.alt, vr.het_or_hom_or_hemi, vr.n_expected_samples))
                    logging.info(", ".join(["%s=%s" % (k, v) for k,v in sorted(counters.items(), key=lambda kv: kv[0])]),)


            else:
                if het_or_hom_or_hemi == "het":
                    if vr.finished and vr.n_expected_samples != min(n_het, MAX_SAMPLES_TO_SHOW_PER_VARIANT):
                        logging.error("het: n expected == %d, n_het in ExAC sites vcf == %d, %s" % (vr.n_expected_samples, n_het, row))
                elif het_or_hom_or_hemi == "hom":
                    if vr.finished and vr.n_expected_samples != min(n_hom, MAX_SAMPLES_TO_SHOW_PER_VARIANT):
                        logging.error("hom: n expected == %d, n_hom in ExAC sites vcf == %d, %s" % (vr.n_expected_samples, n_hom, row))
                elif het_or_hom_or_hemi == "hemi":
                    if vr.finished and vr.n_expected_samples != min(n_hemi, MAX_SAMPLES_TO_SHOW_PER_VARIANT):
                        logging.error("hemi: n expected == %d, n_hemi in ExAC sites vcf == %d, %s" % (vr.n_expected_samples, n_hemi, row))
                else:
                    raise ValueError("Unexpected value for het_or_hom_or_hemi: %s" % str(het_or_hom_or_hemi))
print("Finished")