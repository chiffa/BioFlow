"""
Tests all the functions declared in the utilities module
"""
import os
import unittest
import numpy as np
from itertools import izip

from BioFlow.Utils import Linalg_routines
from BioFlow.Utils import GDF_export
from BioFlow.Utils.GeneralUtils import SanerFilesystem

from BioFlow.Utils.GeneralUtils.SanerFilesystem import wipe_dir, mkdir_recursive


class SanerFilesystemTester(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def just_a_test(self):
        self.assertTrue(True)


class LinalgRoutinesTester(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_sub_matrix_extractor(self):
        self.assertTrue(True)


class GdfExportTester(unittest.TestCase):
    test_location = os.path.join(os.path.dirname(__file__),
                                 'Dumps/GDF_export_test.gdf')
    reference_location = os.path.join(os.path.dirname(__file__),
                                      'UT_examples/GDF_export_reference.gdf')

    @classmethod
    def setUpClass(cls):
        premat = np.zeros((4, 4))

        premat[0, 1] = 1.0
        premat[0, 2] = 4.0
        premat[1, 2] = 0.5
        premat[0, 3] = 0.01

        mkdir_recursive(cls.test_location)

        gdfw = GDF_export.GDF_export_Interface(
            target_fname=cls.test_location,
            field_names=['test'],
            field_types=['VARCHAR'],
            node_properties_dict={'test1': ['test one'],
                                  'test2': ['test two'],
                                  'test3': ['test three'],
                                  'test4': ['test four']},
            mincurrent=0,
            Idx2Label={0: 'test1', 1: 'test2', 2: 'test3', 3: 'test4'},
            Label2Idx={'test1': 0, 'test2': 1, 'test3': 2, 'test4': 3},
            current_Matrix=premat)

        gdfw.write()

    @classmethod
    def tearDownClass(cls):
        wipe_dir('Dumps')

    def test_GDF_export(self):
        with open(self.test_location, 'r') as tested, \
                open(self.reference_location, 'r') as reference:
            for line1, line2 in izip(tested, reference):
                self.assertItemsEqual(line1, line2)

    def test_GDF_exceptions(self):
        pass


if __name__ == "__main__":
    unittest.main()
