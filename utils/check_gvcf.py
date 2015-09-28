import os
import re
import vcf

VARIANT_NOT_IN_ORIGINAL_GVCF = 1  # not sure why this would happen
VARIANT_NOT_IN_NEW_GVCF = 2
VARIANT_GENOTYPES_MISMATCHED = 3


def convert_genotype_to_alleles(vcf_r):
    """Takes a vcf record and converts the genotype string (eg. '0/1') to an
    allele string (eg. 'A/CGG')
    """
    genotype_str = vcf_r.samples[0]["GT"]  # '0/1'
    allele_indexes = map(int, re.split('[/|]', genotype_str))  # [0, 1]
    allele_index_to_genotype_map = [vcf_r.REF] + list(map(str, vcf_r.ALT))  # ['A', 'CGG', 'G']
    genotypes = [allele_index_to_genotype_map[idx] for idx in allele_indexes]  # ['A', 'CGG']
    return "/".join(genotypes)


def check_gvcf(original_gvcf_path, new_gvcf_path, chrom, minrep_pos):
    """Checks whether the GVCF genotype generated by HaplotypeCaller when
    computing the reassemled bam matches the genotype in the original GVCF
    at the given position.

    Don't worry about also checking that they both match the ref/alt allele
    in the final ExAC callset because the representation or allele(s) may have
    changed during the joint calling step.

    Returns: A 3-tuple (success, error_code, error_message) where:
        success = True if the check succeeded
        error_code = An integer error code if the check failed
        error_message = An error message explaining the error code
    """

    #assert os.path.isfile(new_gvcf_path), "new GVCF file not found: %s" % new_gvcf_path

    minrep_pos = int(minrep_pos)


    # being here means a row in the new GVCF matched the chrom/pos of the ExAC call
    variant_found = False
    original_gvcf_file = vcf.Reader(filename=original_gvcf_path)
    for r_original in original_gvcf_file.fetch(chrom, minrep_pos-1, minrep_pos+1):
        if r_original.CHROM == chrom and r_original.POS >= minrep_pos:
            if r_original.POS == minrep_pos:
                variant_found = True
            break
    original_gvcf_file._reader.close()

    if not variant_found:
        return False, VARIANT_NOT_IN_ORIGINAL_GVCF, "%(chrom)s:%(minrep_pos)s not found in original gvcf: %(original_gvcf_path)s" % locals()


    # the new gvcf is small, so just iterate through it until the variant is found
    variant_found = False
    new_gvcf_file = vcf.Reader(filename=new_gvcf_path)
    for r_new in new_gvcf_file:
        if r_new.CHROM == chrom and r_new.POS >= minrep_pos:
            if r_new.POS == minrep_pos:
                variant_found = True
            break
    new_gvcf_file._reader.close()

    if not variant_found:
        return False, VARIANT_NOT_IN_NEW_GVCF, "%(chrom)s:%(minrep_pos)s not found in new gvcf: %(new_gvcf_path)s" % locals()

    # compare genotypes
    original_genotype = convert_genotype_to_alleles(r_original)
    new_genotype = convert_genotype_to_alleles(r_new)
    if original_genotype != new_genotype:
        original_gt = r_original.samples[0]["GT"]
        new_gt = r_new.samples[0]["GT"]
        return False, VARIANT_GENOTYPES_MISMATCHED, "%(chrom)s:%(minrep_pos)s genotypes don't match: %(original_gt)s (%(original_genotype)s) !=  %(new_gt)s (%(new_genotype)s) in   %(original_gvcf_path)s   vs   %(new_gvcf_path)s" % locals()

    return True, 0, None
