'''
Created on Dec 15, 2013
@author: ank
Checks if all the relevant information is in the ConfigParser-defined files
'''
import unittest
from PolyPharma.Utils.ConfigParsers import Configs_common,Configs_parser
# from PolyPharma.Utils.ConfigParsers.Configs_parser import rootdir,configsfiles,shortnames
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

        for ext_DB_type, ext_DB_args in parsedConfs[2].iteritems():
            currentPath = ext_DB_args['location']
            if os.path.isdir(currentPath):
                for fle in [f for f in os.listdir(currentPath) if os.path.isfile(os.path.join(currentPath,f))]:
                    if ext_DB_args['load'] in fle:
                        currentPath=os.path.join(currentPath,fle)
            print ext_DB_type, currentPath
            echo_read(currentPath)



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()