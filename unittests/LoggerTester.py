"""
Tests the logger behavior
"""
import unittest
import os
from BioFlow.configs2 import log_location
from BioFlow.Utils.GeneralUtils.SanerFilesystem import wipe_dir
from BioFlow.Utils.tentative_log_manager import logger


class TestLogs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        logger.debug('This is a debug test')
        logger.info('This is an info logging test')
        logger.warning('This is a warning test')
        logger.error('This is an error test')
        logger.critical('This is a critical test')

    @classmethod
    def tearDownClass(cls):
        wipe_dir(log_location)

    def test_file_creation(self):
        self.assertTrue(os.path.isdir(log_location))

    def test_low_level_logs(self):
        with open(os.path.join(log_location, 'debug.log'), 'r') as source_file:
            contents = ''.join(source_file.read())
            self.assertIn('debug', contents)
            self.assertIn('logging', contents)
            self.assertIn('warning', contents)
            self.assertIn('error', contents)
            self.assertIn('critical', contents)

    def test_high_level_logs(self):
        with open(os.path.join(log_location, 'critical.log'), 'r') as source_file:
            contents = ''.join(source_file.read())
            self.assertIn('critical', contents)
            self.assertNotIn('warning', contents)

    def test_redundant_logs(self):
        with open(os.path.join(log_location, 'info.log'), 'r') as source_file:
            contents = ''.join(source_file.read())
            self.assertEqual(contents.count('debug'), 0)
            self.assertEqual(contents.count('info'), 1)
            self.assertEqual(contents.count('warning'), 1)


if __name__ == "__main__":
    unittest.main()