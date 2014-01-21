'''
Created on Dec 15, 2013
@author: ank
Checks if all the relevant information is in the ConfigParser-defined files
'''
import unittest
import ConfigParser
from PolyPharma.Utils.ConfigParsers.Configs_parser import rootdir,configsfiles,shortnames
import os
from pprint import PrettyPrinter




class ConfigsParser_TestSuite(unittest.TestCase):

    def test_Configs_Parser(self):
        pass
        # print('Testing the contents of the config .ini files')
        # self.assertTrue(os.path.exists(cfgPrs.configsfiles[0]))
        # self.assertTrue(os.path.exists(cfgPrs.configsfiles[0]))
        # pp=PrettyPrinter(indent=4)
        # pp.pprint(cfgPrs.parse_configs())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    generate_configs()
#    unittest.main()