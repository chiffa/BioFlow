"""
Tests all the functions declared in the utilities module
"""
import os
import unittest
import numpy as np
from itertools import izip

from BioFlow.utils import Linalg_routines
from BioFlow.utils import GDF_export
from BioFlow.utils.general_utils import high_level_os_io

from BioFlow.utils.general_utils.high_level_os_io import wipe_dir, mkdir_recursive


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
        cls.test_lapl = Linalg_routines.lil_matrix(np.zeros((4, 4)))
        cls.test_lapl.setdiag([1, 2, 3, 0])
        cls.test_lapl[1, 2] = -2
        cls.test_lapl[2, 1] = -2
        cls.test_lapl[0, 2] = -1
        cls.test_lapl[2, 0] = -1

    @classmethod
    def tearDownClass(cls):
        pass

    def test_normalization(self):
        norm_lapl = Linalg_routines.normalize_laplacian(self.test_lapl).toarray()
        self.assertEqual(norm_lapl[0, 0], 1.)
        self.assertAlmostEqual(norm_lapl[2, 1], -0.816496580)

    def test_eigenvector_analysis(self):
        Linalg_routines.analyze_eigenvects(self. test_lapl, 2, {0: 'zero',
                                                                1: 'one',
                                                                2: 'two',
                                                                3: 'three'})

    def test_clustring(self):
        self.test_lapl[2, 3] = -0.5
        self.test_lapl[3, 2] = -0.5
        idsx = Linalg_routines.cluster_nodes(self.test_lapl).T.tolist()[0]
        self.assertEqual(idsx[1], idsx[2])
        self.assertNotEqual(idsx[0], idsx[3])
        self.assertNotEqual(idsx[1], idsx[0])
        self.assertNotEqual(idsx[1], idsx[3])
        self.assertEqual(Linalg_routines.average_off_diag_in_sub_matrix(self.test_lapl, [1, 2]),
                         -2.)

    def test_average_interset_linkage(self):
        self.assertEqual(Linalg_routines.average_interset_linkage(self.test_lapl,
                                                                  [[0, 1], [2, 3]]), -0.75)


class GdfExportTester(unittest.TestCase):
    test_location = os.path.join(os.path.dirname(__file__),
                                 'dumps/GDF_export_test.gdf')
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
        print cls.test_location
        wipe_dir(cls.test_location)

    def test_GDF_export(self):
        with open(self.test_location, 'r') as tested, \
                open(self.reference_location, 'r') as reference:
            for line1, line2 in izip(tested, reference):
                self.assertItemsEqual(line1, line2)

    def test_GDF_exceptions(self):
        pass


if __name__ == "__main__":
    unittest.main()
