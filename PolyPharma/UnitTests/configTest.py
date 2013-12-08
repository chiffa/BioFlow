'''
Created on Dec 7, 2013

@author: ank
'''
import unittest
from src.config import neo4j


class configTest(unittest.TestCase):

    def test_neo4j_variables_are_declared(self):
        self.assertTrue((neo4j.local_neo4j and neo4j.server_neo4j and neo4j.servers_are_local), "Please declare all the neo4j-related variables")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()