import os
import unittest
os.environ['UNITTESTING'] = 'True'
from unittests.PreProcssingTester import TestRnaCountsProcessor
from unittests.LoggerTester import TestLogs
from unittests.UtilitiesTester import GdfExportTester, LinalgRoutinesTester, SanerFilesystemTester


class HooksConfigTest(unittest.TestCase):

    def test_hooks(self):
        self.assertTrue(True)

    def test_actual_code(self):
        from BioFlow.main_configs import edge_type_filters
        self.assertEqual(edge_type_filters["Group"][0], "is_part_of_collection")


if __name__ == "__main__":
    print TestRnaCountsProcessor.__doc__
    print TestLogs.__doc__
    print GdfExportTester.__doc__
    print LinalgRoutinesTester.__doc__
    print SanerFilesystemTester.__doc__
    unittest.main()
