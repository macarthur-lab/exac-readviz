import os
import vcf
import unittest
from utils.check_gvcf import convert_genotype_to_alleles, check_gvcf



class TestCheckGVCF(unittest.TestCase):

    def setUp(self):
        self.sample1_path = os.path.join(os.path.dirname(__file__), "data/sample1.gvcf.gz")
        self.sample2a_path = os.path.join(os.path.dirname(__file__), "data/sample2a.gvcf.gz")
        self.sample2b_path = os.path.join(os.path.dirname(__file__), "data/sample2b.gvcf.gz")

    def test_convert_genotypes(self):
        # test that all records in sample1.gvcf are converted correctly
        sample1_file = vcf.Reader(filename=self.sample1_path)
        actual = [convert_genotype_to_alleles(r) for r in sample1_file]
        expected = ['G/G', 'C/C', 'T/T', 'T/T', 'C/C', 'T/T', 'T/T', 'A/A', 'T/T',
                    'A/A', 'A/A', 'C/C', 'T/T', 'CTTTTTATTTTTATTTTTA/CTTTTTATTTTTATTTTTA',
                    'A/C', 'TTTTTATTTTTATTTTTA/TTTTTATTTTTATTTTTA', 'A/C', 'T/T', 'A/C',
                    'T/T', 'A/AT', 'T/T', 'T/T', 'A/A', 'C/C', 'T/T', 'A/A', 'A/A', 'G/G',
                    'A/A', 'G/G', 'C/C', 'A/A', 'G/G', 'A/A']

        self.assertListEqual(actual, expected)

        # test that the 4 records in sample2a.gvcf are converted correctly
        sample2a_file = vcf.Reader(filename=self.sample2a_path)
        actual = [convert_genotype_to_alleles(r) for r in sample2a_file]
        expected = ['G/G', 'G/G', 'T/T', 'T/A', 'T/G', 'C/C']
        self.assertListEqual(actual, expected)

        sample2b_file = vcf.Reader(filename=self.sample2b_path)
        actual = [convert_genotype_to_alleles(r) for r in sample2b_file]
        expected = ['G/G', 'G/G', 'G/G', 'G/G', 'T/T', 'G/G', 'G/G', 'G/G',
                    'G/G', 'T/T', 'A/A', 'A/A', 'G/G', 'C/C', 'T/T', 'C/C',
                    'T/T', 'T/A', 'G/G', 'C/G', 'C/G', 'T/T', 'C/C', 'C/C',
                    'C/C', 'C/C', 'G/G', 'G/G', 'A/A', 'T/T', 'G/G', 'G/G',
                    'C/C', 'C/C', 'G/G', 'G/G', 'C/C', 'T/T', 'A/A', 'T/T',
                    'G/G', 'G/G', 'G/G', 'G/G', 'T/T', 'T/T', 'T/T', 'A/A',
                    'T/T', 'G/G', 'G/G', 'C/C', 'G/G', 'G/G', 'C/C', 'G/G',
                    'G/G', 'T/T', 'T/T', 'C/C', 'C/C', 'A/A', 'C/C', 'G/G',
                    'C/C', 'C/C', 'T/T', 'T/T', 'G/G']

    def test_check_gvcf(self):
        # sanity check
        for path in [self.sample1_path, self.sample2a_path, self.sample2b_path]:
            for r in vcf.Reader(filename=path):
                success, error_code, error_message = check_gvcf(path, path, r.CHROM, r.POS)
                self.assertTrue(success, "check_gvcf against self failed at %(r)s in %(path)s: %(error_code)s %(error_message)s " % locals())

        # check for matches
        for r in vcf.Reader(filename=self.sample2b_path):
            if r.CHROM != "2" or r.POS < 905393 or r.POS > 905494:
                continue

            success, error_code, error_message = check_gvcf(self.sample2a_path, self.sample2b_path, r.CHROM, r.POS)

            #print(r.CHROM, r.POS, success, error_message)

            if r.POS == 905394: self.assertTrue(success)
            if r.POS == 905492: self.assertTrue(success)
            if r.POS == 905493: self.assertFalse(success)
