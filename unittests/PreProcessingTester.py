"""
Tests the pre-processing module
"""
import os
os.environ['UNITTESTING'] = 'True'
import unittest
import numpy as np

from bioflow.pre_processing import rna_counts_analysis as RCP


class TestRnaCountsProcessor(unittest.TestCase):
    counts_ref = os.path.join(os.path.dirname(__file__),
                              'UT_examples/counts.tsv')

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_table_loading(self):
        names, lengths, data = RCP.load_rna_counts_table(self.counts_ref, 6)
        self.assertEqual(data.shape, (12, 6))
        self.assertEqual(names[0, 0], 'ENSMUSG00000000001')
        self.assertEqual(names[-1, 0], 'ENSMUSG00000000093')
        self.assertEqual(lengths[1], 902)

    def test_counts_filter(self):
        _, _, table = RCP.load_rna_counts_table(self.counts_ref, 6)
        comp_array = RCP.counts_filter(table, [[0, 1, 2], [3, 4, 5]], 5)
        ref_array = np.array([1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]).astype(np.bool_)
        self.assertCountEqual(comp_array, ref_array)

    def test_rpkm_conversion(self):
        _, lengths, data = RCP.load_rna_counts_table(self.counts_ref, 6)
        comp_data = RCP.convert_to_rpkm(lengths, data)
        self.assertAlmostEqual(comp_data[0, 0]/1e+8, 3.07044661)
        self.assertAlmostEqual(comp_data[-1, -1]/1e+8, 2.79674060)

    def test_stats(self):
        names, lengths, data = RCP.load_rna_counts_table(self.counts_ref, 6)
        filter_mask = RCP.counts_filter(data, [[0, 1, 2], [3, 4, 5]], 5)
        data = data[filter_mask, :]
        names, lengths = (names[filter_mask, :], lengths[filter_mask, :])
        rpkms = RCP.convert_to_rpkm(lengths, data)
        comp_val = RCP.significantly_different_genes(rpkms,
                                                     [[0, 1, 2], [3, 4, 5]],
                                                     [[0, 1]], 0.3)[0][1]
        ref_val = np.array([1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]).astype(np.bool_)
        self.assertCountEqual(comp_val, ref_val)

    def test_suite(self):
        exp_groups = [[0, 1, 2], [3, 4, 5]]
        test_groups_to_compare = [[0, 1]]
        comp_val = RCP.run_analysis_suite(self.counts_ref, 6, exp_groups,
                                     test_groups_to_compare, 5, 0.3)[0]
        ref_val = np.array([1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]).astype(np.bool_)
        self.assertCountEqual(comp_val, ref_val)

if __name__ == "__main__":
    unittest.main()
