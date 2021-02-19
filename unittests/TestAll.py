import os
import unittest
os.environ['UNITTESTING'] = 'True'
from unittests.PreProcessingTester import TestRnaCountsProcessor
from unittests.LoggerTester import TestLogs
from unittests.UtilitiesTester import GdfExportTester, LinalgRoutinesTester, SanerFilesystemTester
from unittests.ParserTester import GoParserTester, UniprotParserTester, ReactomeParseTester
from unittests.ConductionTester import ConductionRoutinesTester


class HooksConfigTest(unittest.TestCase):

    def test_hooks(self):
        self.assertTrue(True)


if __name__ == "__main__":
    my_list = [
        TestRnaCountsProcessor.__doc__, TestLogs.__doc__, GdfExportTester.__doc__,
        LinalgRoutinesTester.__doc__, SanerFilesystemTester.__doc__, GoParserTester.__doc__,
        UniprotParserTester.__doc__, ReactomeParseTester.__doc__,
        ConductionRoutinesTester.__doc__]
    unittest.main()
