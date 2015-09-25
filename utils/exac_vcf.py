"""
Utility methods for parsing the full ExAC vcf.

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

import logging
import re

def create_vcf_row_parser(header_line, valid_sample_ids):
    """Defines and returns a function that can parse a single VCF row
    Args:
      header_line: The last line of the VCF header - the one that defines columns.
      valid_sample_ids: A set of valid sample ids.
    """
    header_fields = header_line.strip("\n").split("\t")

    assert header_fields[0] == "#CHROM" and header_fields[1] == "POS", \
        "Unexpected header_fields 1: %s" % str(header_fields[0:9])
    assert header_fields[3] == "REF" and header_fields[4] == "ALT", \
        "Unexpected header_fields 2: %s" % str(header_fields[0:9])
    assert header_fields[8] == "FORMAT", \
        "Unexpected header_fields 3: %s" % str(header_fields[0:9])

    # sanity check for sample_ids
    sample_ids = header_fields[9:]
    for sample_id in sample_ids:
        if sample_id not in valid_sample_ids:
            logging.error("ERROR: vcf sample id '%s' is not in the vcf info table" % sample_id)

    expected_genotype_format = "GT:AD:DP:GQ:PL"
    GT_idx = expected_genotype_format.split(":").index("GT")
    GQ_idx = expected_genotype_format.split(":").index("GQ")
    DP_idx = expected_genotype_format.split(":").index("DP")

    def vcf_row_parser(fields):
        """Takes a tuple of column values from a VCF row and return a tuple of (chrom, pos, ref, alt_alleles, n_het, n_hom, genotypes)
        where alt_alleles, n_het, and n_hom are lists of strings, and genotypes is a dictionary that maps sample_id to a 4-tuple:
        (gt_ref, gt_alt, GQ, DP)
        """
        chrom, pos, ref, alt_list, info = fields[0], fields[1], fields[3], fields[4].split(","), fields[7]

        m = re.search("AC_Adj=([\d,]+);.*AC_Hom=([\d,]+);", info)
        assert m, "%s:%s: AC_Adj=([\d,]+);AC_Hom=([\d,]+); not found in %s" % (chrom, pos, info)

        ac_adj_list = list(map(int, m.group(1).split(",")))
        n_hom_list = list(map(int, m.group(2).split(",")))
        n_het_list = list(map(lambda t: t[0] - t[1]*2, zip(ac_adj_list, n_hom_list)))

        assert len(n_het_list) == len(alt_list), "%s:%s: Unexpected number of AC_Het: %s - %s" % (chrom, pos, n_het_list, "\t".join(fields[0:9]))
        assert len(n_hom_list) == len(alt_list), "%s:%s: Unexpected number of AC_Hom: %s - %s" % (chrom, pos, n_hom_list, "\t".join(fields[0:9]))

        assert fields[8] == expected_genotype_format, \
            "%s:%s: Unexpected genotype format: '%s'. Expected: %s in %s" % (chrom, pos, fields[8], expected_genotype_format, fields[0:8])

        genotypes = fields[9:]
        assert len(sample_ids) == len(genotypes), \
            "%s:%s: Unexpected num sample ids (%s) vs num genotypes (%s) in %s" % (chrom, pos, len(sample_ids), len(genotypes), fields[0:8])

        sample_id_to_genotype = {}
        for sample_id, genotype in zip(sample_ids, genotypes):
            genotype_values = genotype.split(":")
            GT = genotype_values[GT_idx]
            if GT == "./.":
                gt_ref = gt_alt = None
                GQ = DP = None
            else:
                gt_ref, gt_alt = GT.split("/")
                try:
                    gt_ref = int(gt_ref)
                    gt_alt = int(gt_alt)

                except ValueError:
                    logging.error("ERROR: %s:%s - couldn't parse genotype %s in %s" % (chrom, pos, GT, fields[0:8]))

                assert gt_ref <= len(alt_list) and gt_alt <= len(alt_list), "ERROR: %s:%s - genotype numbers %s out of bounds in %s" % (chrom, pos, GT, fields[0:8])

                try:
                    GQ = float(genotype_values[GQ_idx])
                    DP = float(genotype_values[DP_idx])
                except ValueError:
                    logging.error("ERROR: %s:%s - couldn't parse %s genotype: %s in %s" % (chrom, pos, sample_id, genotype, fields[0:8]))
                    # This error happens for genotypes like 1/1:0,0:.:6:70,6,0.
                    # set GQ, DP = 0 so this sample will be filtered out by the GQ<20,DP<10 filter
                    GQ = DP = 0

            sample_id_to_genotype[sample_id] = (gt_ref, gt_alt, GQ, DP)


        return chrom, pos, ref, alt_list, n_het_list, n_hom_list, sample_id_to_genotype

    return vcf_row_parser

