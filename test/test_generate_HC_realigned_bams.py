import unittest
import generate_HC_realigned_bams

class ChooseSamplesTest(unittest.TestCase):

    def setUp(self):

        self.genotypes = {
            #  sample name : (gt_type, GQ)
                    "s1het": (1, 10),
                    "s2het": (1, 20),
                    "s3het": (1, 30),
                    "s4het": (1, 40),
                    "s5het": (1, 50),
                    "s6het": (1, 60),

                    "s1hom": (2, 10),
                    "s2hom": (2, 20),
                    "s3hom": (2, 30),
                    "s4hom": (2, 40),
                    "s5hom": (2, 50),
                    "s6hom": (2, 60),
        }

        # stub out the get_bam_path method to only test one method, and get rid of checking whether the bam exists on disk
        def get_bam_path_method_stub(sample_name, sample_name_to_bam_path):
            return sample_name+".bam"

        self.original_get_bam_paths_method = generate_HC_realigned_bams.get_bam_path
        generate_HC_realigned_bams.get_bam_path = get_bam_path_method_stub

        self.sample_name_include_status = {}
        for sample_name in self.genotypes:
            self.sample_name_include_status[sample_name] = True

        self.sample_name_to_bam_path = {}
        self.all_sample_bam_paths = []
        for sample_name in self.genotypes:
            bam_path = sample_name+".bam"
            self.sample_name_to_bam_path[sample_name] = bam_path
            self.all_sample_bam_paths.append(bam_path)


    def test_choose_samples(self):
        # test het
        bam_paths = generate_HC_realigned_bams.choose_samples("het", self.genotypes, self.sample_name_include_status, self.sample_name_to_bam_path)
        expected_bam_paths = ['s6het.bam', 's5het.bam', 's4het.bam', 's3het.bam', 's2het.bam']
        self.assertListEqual(bam_paths, expected_bam_paths)

        # test hom
        bam_paths = generate_HC_realigned_bams.choose_samples("hom", self.genotypes, self.sample_name_include_status, self.sample_name_to_bam_path)
        expected_bam_paths = ['s6hom.bam', 's5hom.bam', 's4hom.bam', 's3hom.bam', 's2hom.bam']
        self.assertListEqual(bam_paths, expected_bam_paths)

        # delete every other hom genotype so there are < 5 homs
        del self.genotypes["s2hom"]
        del self.genotypes["s4hom"]
        del self.genotypes["s6hom"]
        bam_paths = generate_HC_realigned_bams.choose_samples("hom", self.genotypes, self.sample_name_include_status, self.sample_name_to_bam_path)
        expected_bam_paths = ['s5hom.bam', 's3hom.bam', 's1hom.bam']
        self.assertListEqual(bam_paths, expected_bam_paths)

        # delete every other het genotype so there are < 5 hets
        del self.genotypes["s2het"]
        del self.genotypes["s4het"]
        del self.genotypes["s6het"]
        bam_paths = generate_HC_realigned_bams.choose_samples("het", self.genotypes, self.sample_name_include_status, self.sample_name_to_bam_path)
        expected_bam_paths = ['s5het.bam', 's3het.bam', 's1het.bam']
        self.assertListEqual(bam_paths, expected_bam_paths)

        # delete remaining hets
        del self.genotypes["s1hom"]
        del self.genotypes["s3hom"]
        del self.genotypes["s5hom"]
        bam_paths = generate_HC_realigned_bams.choose_samples("hom", self.genotypes, self.sample_name_include_status, self.sample_name_to_bam_path)
        self.assertListEqual(bam_paths, [])

    def tearDown(self):
        generate_HC_realigned_bams.get_bam_path = self.original_get_bam_paths_method
