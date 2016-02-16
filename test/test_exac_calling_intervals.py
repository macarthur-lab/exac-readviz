import os
import peewee
import unittest
from playhouse.test_utils import test_database
from utils.exac_calling_intervals import *

test_db = peewee.SqliteDatabase(':memory:')

class TestExacCallingIntervals(unittest.TestCase):

    def load_test_data(self):
        parse_exac_calling_intervals(
            os.path.join(os.path.dirname(__file__),
                         "data/calling_regions.interval_list"))

    def test_exac_calling_intervals(self):
        with test_database(test_db, (ExacCallingInterval,)):
            self.load_test_data()   # this data will be loaded into `test_db`

            # test that data was loaded correctly
            intervals = list(ExacCallingInterval.select())
            self.assertEqual(str(intervals[0]), "1:12141-12277")
            self.assertEqual(str(intervals[1]), "1:12546-12771")
            self.assertEqual(str(intervals[-1]), "MT:12287-15937")

            # test get_overlapping_calling_interval
            intervals = list(ExacCallingInterval.select())
            self.assertEqual(str(get_overlapping_calling_interval("1", 12546)), "1:12546-12771")
            self.assertEqual(str(get_overlapping_calling_interval("1", 12671)), "1:12546-12771")
            self.assertEqual(str(get_overlapping_calling_interval("1", 12771)), "1:12546-12771")
            self.assertRaises(ValueError, get_overlapping_calling_interval, "1", 12545)
            self.assertRaises(ValueError, get_overlapping_calling_interval, "1", 12772)

            # test get_adjacent_calling_interval
            l,i,r = get_adjacent_calling_intervals("1", 12546)
            self.assertEqual(len(l), 1)
            self.assertEqual(len(r), 1)
            self.assertListEqual(list(map(str, [l[0],i,r[0]])),
                [ "1:12141-12277", "1:12546-12771", "1:13354-13689"])

            l,i,r = get_adjacent_calling_intervals("1", 12546, n_left=2, n_right=2)
            self.assertEqual(len(l), 1)
            self.assertEqual(len(r), 2)
            self.assertListEqual(list(map(str, [l[0], i, r[0], r[1]])),
                [ "1:12141-12277", "1:12546-12771", "1:13354-13689", "1:17319-17486"])
