'''
Created on Dec 15, 2013
@author: ank
Checks if all the relevant information is in the ConfigParser-defined files
'''
import unittest
from PolyPharma.Utils.ConfigParsers import Configs_common,Configs_parser
import os
from pprint import PrettyPrinter

def echo_read(path_to_file):
    fle=open(path_to_file,'r')
    for line in fle.readlines()[:10]:
        print line
    fle.close()

class ConfigsParser_TestSuite(unittest.TestCase):

    def test_Configs_Parser(self):
        print('Testing the contents of the config .ini files')
        for i in range(0,3):
            self.assertTrue(os.path.exists(Configs_common.configsfiles[i]))
        parsedConfs = Configs_parser.parse_configs()

        pp = PrettyPrinter(indent = 4)
        pp.pprint(parsedConfs)

        ReadFiles=Configs_parser.sourcefile_compilator(parsedConfs[2])
        pp.pprint(ReadFiles)
        for path in ReadFiles.values():
            echo_read(path)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()