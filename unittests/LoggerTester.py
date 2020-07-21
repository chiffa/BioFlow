"""
Tests the logger behavior
"""
import unittest
import os
from bioflow.main_configs import log_location
from bioflow.utils.general_utils.high_level_os_io import wipe_dir
from bioflow.utils.log_behavior import get_logger

log = get_logger('LoggerTester @ unittests')


class TestLogs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        log.debug('This is a debug test')
        log.info('This is an info logging test')
        log.warning('This is a warning test')
        log.error('This is an error test')
        log.critical('This is a critical test')

    @classmethod
    def tearDownClass(cls):
        try:
            wipe_dir(log_location)
        except OSError:
            pass

    def test_file_creation(self):
        self.assertTrue(os.path.isdir(log_location))

    def test_low_level_logs(self):
        with open(os.path.join(log_location, 'debug.log'), 'rt') as source_file:
            contents = ''.join(source_file.read())
            self.assertIn('debug', contents)
            self.assertIn('logging', contents)
            self.assertIn('warning', contents)
            self.assertIn('error', contents)
            self.assertIn('critical', contents)

    def test_high_level_logs(self):
        with open(os.path.join(log_location, 'critical.log'), 'rt') as source_file:
            contents = ''.join(source_file.read())
            self.assertIn('critical', contents)
            self.assertNotIn('warning', contents)

    def test_redundant_logs(self):
        with open(os.path.join(log_location, 'info.log'), 'rt') as source_file:
            contents = ''.join(source_file.read())
            self.assertEqual(contents.count('debug'), 0)
            self.assertEqual(contents.count('info'), 1)
            self.assertEqual(contents.count('warning'), 1)


if __name__ == "__main__":
    unittest.main()
