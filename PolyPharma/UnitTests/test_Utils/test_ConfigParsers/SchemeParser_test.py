'''
Created on Dec 26, 2013

@author: ank
'''
import unittest


class Test(unittest.TestCase):
# TODO: check for validation and failure in case of every case in the 
# inner functions of the large functions

    def setUp(self):
        # a single pass test;
        # multiple fail tests
            # - IO fails
            # - altered core classes
            # - lacking inheritance
            # - additional reference to a baseclass
            # - circular reference
            # - wrong type declaration
        pass


    def tearDown(self):
        # delete all the IO files created by the setUp module
        pass


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()