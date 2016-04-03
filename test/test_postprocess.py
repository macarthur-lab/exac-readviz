import unittest
from utils.postprocess_reassembled_bam import do_intervals_intersect, interval_union


class TestPostProcess(unittest.TestCase):

    def test_interval_unionc(self):
        self.assertEqual(interval_union((1,5), (3,4)),    (1,5))
        self.assertEqual(interval_union((1,5), (1e9, 0)), (1,5))
        self.assertEqual(interval_union((1,5), (5, 10)),  (1,10))

    def test_interval_intersection(self):
        self.assertTrue(do_intervals_intersect((1,5), (2,8)))

        self.assertFalse(do_intervals_intersect((1,5), (5,8)))
        self.assertFalse(do_intervals_intersect((5,8), (1,5)))

        self.assertTrue(do_intervals_intersect((6,8), (1,10)))
