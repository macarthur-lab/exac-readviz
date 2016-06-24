"""
Utility methods that can be used for parsing the full ExAC genotypes vcf or the sites VCF.

The VCF must have the following INFO fields:

AC_Adj=([\d,]+);
AC_Hom=([\d,]+);
AC_Hemi=([\d,]+);   (on chrom X, Y)

The VCF genotype format must be:

GT:AD:DP:GQ:PL


Example usage:

   # open full ExAC VCF
   tabix_file = pysam.TabixFile(filename=args.full_vcf, parser=pysam.asTuple())

   # use VCF header line to create row parser function
   last_header_line = list(tabix_file.header)[-1].decode("utf-8", "ignore")
   parse_vcf_row = create_vcf_row_parser(last_header_line, set(sample_id_include_status.keys()))

   for row in vcf_iterator:
        chrom, pos, ref, alt_alleles, genotypes = parse_vcf_row(row)

        # iterate over alt alleles (in case this row is multi-allelic)
        for alt_allele_index, alt in enumerate(alt_alleles):
            # minrep
            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(pos, ref, alt)

"""

import collections
import logging
import re


from utils.minimal_representation import get_minimal_representation


EXPECTED_GENOTYPE_FORMAT = "GT:AD:DP:GQ:PL"
GENOTYPE_FORMAT_FIELDS = EXPECTED_GENOTYPE_FORMAT.split(":")

GT_idx = GENOTYPE_FORMAT_FIELDS.index("GT")
AD_idx = GENOTYPE_FORMAT_FIELDS.index("AD")
DP_idx = GENOTYPE_FORMAT_FIELDS.index("DP")
GQ_idx = GENOTYPE_FORMAT_FIELDS.index("GQ")


def parse_info_field(fields, alt_alleles, info_field, chrom = None, pos = None):
    """Parses a VCF info fields and returns a tuple with 3 lists: n_het_list, n_hom_list, n_hemi_list"""

    m = re.search("AC_Adj=([\d,]+);.*AC_Hom=([\d,]+);", info_field)
    assert m, "%s:%s: AC_Adj=([\d,]+);AC_Hom=([\d,]+); not found in %s" % (chrom, pos, info_field)

    ac_adj_list = tuple(map(int, m.group(1).split(",")))

    # list with one element per allele. Each element is the number of homozygous samples expected for that allele
    n_hom_list = tuple(map(int, m.group(2).split(",")))

    # list with one element per allele. Each element is the number of hemizygous samples expected for that allele
    if chrom in ("X", "Y"):
        m = re.search("AC_Hemi=([\d,]+);", info_field)
        assert m, "%s:%s: AC_Hemi=([\d,]+); not found in %s" % (chrom, pos, info_field)

        n_hemi_list = tuple(map(int, m.group(1).split(",")))
    else:
        n_hemi_list = tuple([0]*len(ac_adj_list))

    # list with one element per allele. Each element is the number of het samples expected for that allele
    n_het_list = list(map(lambda t: t[0] - t[1]*2 - t[2], zip(ac_adj_list, n_hom_list, n_hemi_list)))

    assert len(n_het_list) == len(alt_alleles), "%s:%s: Unexpected number of AC_Het: %s - %s" % (chrom, pos, n_het_list, "\t".join(fields[0:9]) if fields else "")
    assert len(n_hom_list) == len(alt_alleles), "%s:%s: Unexpected number of AC_Hom: %s - %s" % (chrom, pos, n_hom_list, "\t".join(fields[0:9]) if fields else "")
    assert len(n_hemi_list) == len(alt_alleles), "%s:%s: Unexpected number of AC_Hemi: %s - %s" % (chrom, pos, n_hemi_list, "\t".join(fields[0:9]) if fields else "")

    return n_het_list, n_hom_list, n_hemi_list


def parse_genotypes(fields, alt_alleles, vcf_sample_ids, chrom = None, pos = None):
    """Parses all genotypes and returns a dictionary that maps each sample id to its Genotype"""

    # define Genotype tuple
    Genotype = collections.namedtuple("Genotype", [
        "gt_ref",
        "gt_alt",
        "AD",
        "DP",
        "GQ"
    ])

    assert fields[8] == EXPECTED_GENOTYPE_FORMAT, \
        "%s:%s: Unexpected genotype format: '%s'. Expected: %s in %s" % (chrom, pos, fields[8], EXPECTED_GENOTYPE_FORMAT, fields[0:8])

    genotypes = fields[9:]
    assert len(vcf_sample_ids) == len(genotypes), \
        "%s:%s: Unexpected num sample ids (%s) vs num genotypes (%s) in %s" % (chrom, pos, len(vcf_sample_ids), len(genotypes), fields[0:8])

    sample_id_to_genotype = {}
    for vcf_sample_id, genotype in zip(vcf_sample_ids, genotypes):
        genotype_values = genotype.split(":")
        GT = genotype_values[GT_idx]
        if GT == "./.":
            gt_ref = gt_alt = None
            AD = DP = GQ = None
        else:
            gt_ref, gt_alt = GT.split("/")
            try:
                gt_ref = int(gt_ref)
                gt_alt = int(gt_alt)

                assert gt_ref <= len(alt_alleles) and gt_alt <= len(alt_alleles), "ERROR: %s:%s - genotype numbers %s out of bounds in %s" % (chrom, pos, GT, fields[0:8])

            except ValueError:
                logging.error("ERROR: %s:%s - couldn't parse genotype %s in %s" % (chrom, pos, GT, fields[0:8]))
                gt_ref = gt_alt = None

            try:
                AD = list(map(int, genotype_values[AD_idx].split(",")))
                DP = float(genotype_values[DP_idx])
                GQ = float(genotype_values[GQ_idx])
            except ValueError:
                logging.error("ERROR: %s:%s - couldn't parse %s genotype: %s in %s" % (chrom, pos, vcf_sample_id, genotype, fields[0:8]))
                # This error happens for genotypes like 1/1:0,0:.:6:70,6,0.
                # set GQ, DP = 0 so this sample will be filtered out by the GQ<20,DP<10 filter
                AD = None
                DP = GQ = 0

        sample_id_to_genotype[vcf_sample_id] = Genotype(gt_ref, gt_alt, AD, DP, GQ)

    return sample_id_to_genotype


def create_vcf_row_parser(header_line, valid_sample_ids=None):
    """Defines and returns a function that can parse a single VCF row (represented by a tuple of column values) and
    yields one or more VariantGT objects as defined below.

    Args:
        header_line: The last line of the VCF header - the one that defines columns.
        valid_sample_ids: A set of valid sample ids. If specified, it should contain all sample ids that appear in the VCF,
            of should be None in which case all genotypes will be ignored.
    Return:
        Function that parses vcf row and yields one or more VariantGT objects.
    """

    # define object returned by vcf_row_to_variants
    VariantGT = collections.namedtuple('VariantGT', [
        'chrom',   # chromosome (eg. '1', 'X', etc.)
        'pos',     # minreped pos
        'ref',     # minreped ref allele
        'alt',     # minreped pos allele
        'alt_allele_index',     # index of this alt allele (eg. 0 for the 1st alt allele, etc.)
        'het_or_hom_or_hemi',   # string value - one of:  "het", "hom", "hemi"
        'n_expected_samples',   # number of samples expected to be het_or_hom_or_hemi based on AC_
        'all_genotypes_in_row'  # dictionary that maps sample_id to a 5-tuple:(gt_ref, gt_alt, AD, DP, GQ)
    ])


    header_fields = header_line.strip("\n").split("\t")

    assert header_fields[0] == "#CHROM" and header_fields[1] == "POS", \
        "Unexpected header_fields 1: %s" % str(header_fields[0:9])
    assert header_fields[3] == "REF" and header_fields[4] == "ALT", \
        "Unexpected header_fields 2: %s" % str(header_fields[0:9])
    if valid_sample_ids is not None:
        assert header_fields[8] == "FORMAT", \
            "Unexpected header_fields 3: %s" % str(header_fields[0:9])

    # sanity check for sample_ids
    if valid_sample_ids is not None:
        vcf_sample_ids = header_fields[9:]
        for vcf_sample_id in vcf_sample_ids:
            if vcf_sample_id not in valid_sample_ids:
                logging.error("ERROR: vcf sample id '%s' is not in the set of %d valid_sample_ids" % (
                    vcf_sample_id, len(valid_sample_ids)))


    def vcf_row_to_variants(fields, het_or_hom_or_hemi=None):
        """Takes a single VCF row (represented by a tuple of column values) and yields one or more VariantGT objects as defined above.

        Args:
            fields: tuple of VCF column values
            het_or_hom_or_hemi: (optional) "het" or "hom" or "hemi" string. If not None, only this type of variant will be yielded.
        """

        chrom, pos, ref, alt_alleles, info_field = fields[0], fields[1], fields[3], fields[4].split(","), fields[7]

        if het_or_hom_or_hemi is not None:
            possible_genotypes = (het_or_hom_or_hemi, )
        else:
            possible_genotypes = ("het", "hom", "hemi") if chrom in ('X', 'Y') else ("het", "hom")

        # parse info field
        n_het_list, n_hom_list, n_hemi_list = parse_info_field(fields, alt_alleles, info_field, chrom=chrom, pos=pos)

        # handle genotypes if valid_sample_ids arg is not None
        sample_id_to_genotype = None
        if valid_sample_ids is not None:
            sample_id_to_genotype = parse_genotypes(fields, alt_alleles, vcf_sample_ids, chrom=chrom, pos=pos)


        for alt_allele_index, (alt, n_het, n_hom, n_hemi) in enumerate(zip(alt_alleles, n_het_list, n_hom_list, n_hemi_list)):
            minrep_pos, minrep_ref, minrep_alt = get_minimal_representation(pos, ref, alt)

            for het_or_hom_or_hemi in possible_genotypes:
                if het_or_hom_or_hemi == "het":
                    n_expected_samples = n_het
                elif het_or_hom_or_hemi == "hom":
                    n_expected_samples = n_hom
                elif het_or_hom_or_hemi == "hemi":
                    n_expected_samples = n_hemi
                else:
                    raise ValueError("Unexpected value for het_or_hom_or_hemi: %s" % str(het_or_hom_or_hemi))

                yield VariantGT(chrom=chrom,
                                pos=minrep_pos,
                                ref=minrep_ref,
                                alt=minrep_alt,
                                alt_allele_index=alt_allele_index,
                                het_or_hom_or_hemi=het_or_hom_or_hemi,
                                n_expected_samples=n_expected_samples,
                                all_genotypes_in_row=sample_id_to_genotype)

    return vcf_row_to_variants


def create_variant_iterator_from_vcf(vcf_row_iterator, parse_vcf_row):
    """Iterate over & parse the vcf rows, yielding one or more VariantGT objects from each row"""

    for row in vcf_row_iterator:
        for variant in parse_vcf_row(row):
            yield variant
