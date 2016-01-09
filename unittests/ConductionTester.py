import os
import unittest
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
import warnings
from bioflow.algorithms_bank import conduction_routines as cr


class ConductionRoutinesTester(unittest.TestCase):

    ref_results = os.path.join(os.path.dirname(__file__), 'UT_examples/conduction_reference.dmp')

    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Changing the sparsity structure")
            cls.test_laplacian = csc_matrix(np.zeros((4, 4)))
            cls.test_laplacian.setdiag([1, 2, 3, 0])
            cls.test_laplacian[1, 2] = -2
            cls.test_laplacian[2, 1] = -2
            cls.test_laplacian[0, 2] = -1
            cls.test_laplacian[2, 0] = -1

    def test_sparse_abs(self):
        ref = np.abs(self.test_laplacian.toarray())
        calc = cr.sparse_abs(self.test_laplacian).toarray()
        self.assertListEqual(calc.T.tolist(), ref.T.tolist())

    def test_build_sink_source_current_array(self):
        ref = np.array([[1, -1, 0, 0]]).T
        calc = cr.build_sink_source_current_array((0, 1), (4, 4))
        self.assertListEqual(calc.T.tolist(), ref.T.tolist())

    def test_build_potentials(self):
        calc = cr.get_potentials(self.test_laplacian, (0, 3))
        self.assertTrue(calc[0, 0] > 1e5)
        calc = cr.get_potentials(self.test_laplacian, (0, 1))
        ref = np.array([[0.833333333, -0.666666666, -0.166666666666, 0.0]]).T
        self.assertTrue(np.mean(np.abs(calc - ref)) < 1e-9)

    def test_get_current_matrix(self):
        potentials = cr.get_potentials(self.test_laplacian, (0, 1))
        currents, triu_currents = cr.get_current_matrix(self.test_laplacian, potentials)
        calc = triu_currents.toarray()[:, 2].tolist()
        ref = np.array([-1, 1, 0, 0])
        self.assertTrue(np.mean(np.abs(calc - ref)) < 1e-9)

    def test_get_current_through_nodes(self):
        potentials = cr.get_potentials(self.test_laplacian, (0, 1))
        currents, triu_currents = cr.get_current_matrix(self.test_laplacian, potentials)
        triu_currents = lil_matrix(triu_currents)
        node_currents = cr.get_current_through_nodes(triu_currents)
        ref = np.array([1, 1, 1, 0])
        self.assertTrue(np.mean(np.abs(node_currents - ref)) < 1e-9)

    def test_group_induced_current(self):
        triu_currents = cr.group_edge_current(self.test_laplacian, [0, 1, 2])
        calc = triu_currents.toarray()[:, 2].tolist()
        ref = np.array([1.6666666666, 2.666666666, 0, 0])
        self.assertTrue(np.mean(np.abs(calc - ref)) < 1e-9)

    def test_group_induced_current_memoized(self):
        triu_currents, memoizer = cr.group_edge_current_memoized(self.test_laplacian, [0, 1, 2])
        calc = triu_currents.toarray()[:, 2].tolist()
        ref = np.array([1.6666666666, 2.666666666, 0, 0])/3.
        self.assertTrue(np.mean(np.abs(calc - ref)) < 1e-9)

        chm1 = np.zeros((4, 4))
        chm2 = np.zeros((4, 4))
        chm3 = np.zeros((4, 4))
        chm1[0, 2] = -1
        chm1[1, 2] = 1
        chm2[1, 2] = -1
        chm3[0, 2] = -1

        calc = memoizer[(0, 1)][1].toarray()
        self.assertTrue(np.mean(np.abs(calc - chm1)) < 1e-9)
        calc = memoizer[(1, 2)][1].toarray()
        self.assertTrue(np.mean(np.abs(calc - chm2)) < 1e-9)
        calc = memoizer[(0, 2)][1].toarray()
        self.assertTrue(np.mean(np.abs(calc - chm3)) < 1e-9)

    def test_laplacian_reachable_filter(self):
        chm = np.zeros((4, 4))
        chm[0, 0] = 1
        chm[2, 2] = 1
        chm[2, 0] = -1
        chm[0, 2] = -1

        calc_m = cr.laplacian_reachable_filter(self.test_laplacian, (0, 2)).toarray()
        self.assertTrue(np.mean(np.abs(calc_m - chm)) < 1e-9)

    def test_current_with_reach_limitations(self):
        calc_m, current = cr.group_edge_current_with_limitations(self.test_laplacian,
                                                                 (0, 2), [0, 2])
        chm = np.zeros((4, 4))
        chm[0, 2] = -1
        calc_m = calc_m.toarray()
        self.assertTrue(np.mean(np.abs(calc_m - chm)) < 1e-9)


if __name__ == "__main__":
    unittest.main()
