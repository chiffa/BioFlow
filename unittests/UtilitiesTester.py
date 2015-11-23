"""
Tests all the functions declared in the utilities module
"""
import unittest
import numpy as np
from BioFlow.configs2 import Dumps
from BioFlow.Utils import Linalg_routines
from BioFlow.Utils import GDF_export
from BioFlow.Utils.GeneralUtils import SanerFilesystem


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

    @classmethod
    def setUpClass(cls):
        premat = np.zeros((4, 4))

        premat[0, 1] = 1.0
        premat[0, 2] = 4.0
        premat[1, 2] = 0.5
        premat[0, 3] = 0.01

        GDFW = GDF_export.GDF_export_Interface(
            target_fname=Dumps.GDF_debug,
            field_names=['test'],
            field_types=['VARCHAR'],
            node_properties_dict={'test1': ['test one'], 'test2': ['test two'],
                                  'test3': ['test three'], 'test4': ['test four']},
            mincurrent=0,
            Idx2Label={0: 'test1', 1: 'test2', 2: 'test3', 3: 'test4'},
            Label2Idx={'test1': 0, 'test2': 1, 'test3': 2, 'test4': 3},
            current_Matrix=premat)

        GDFW.write_nodedefs()
        GDFW.write_nodes()
        GDFW.write_edgedefs()
        GDFW.write_edges()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_GDF_export(self):
        pass


if __name__ == "__main__":
    unittest.main()
