import os
import unittest
from BioFlow.bio_db_parsers.gene_ontology_parser import GOTermsParser


class GoParserTester(unittest.TestCase):

    ref_obo = os.path.join(os.path.dirname(__file__), 'UT_examples/test_go.obo')

    @classmethod
    def setUpClass(cls):
        cls.terms, cls.term_rels = GOTermsParser().parse_go_terms(cls.ref_obo)

    def test_proper_parsing(self):
        self.assertIn('0000001', self.terms.keys())
        self.assertIn('0000002', self.terms.keys())

    def test_namespaces(self):
        self.assertEqual('biological_process', self.terms['0000001']['namespace'])

    def test_names(self):
        self.assertEqual('mitochondrion inheritance', self.terms['0000001']['name'])

    def test_relations(self):
        self.assertIn(('0000001', 'is_a', '0048311'), self.term_rels)
        self.assertIn(('0000001', 'is_a', '0048308'), self.term_rels)
        self.assertIn(('0000002', 'is_a', '0007005'), self.term_rels)

    def test_obsolescence(self):
        self.assertNotIn('0000039', self.terms.keys())

    def test_chebi_parsing(self):
        self.assertIn('CHEBI', self.terms['0000036'].keys())
        self.assertIn('22221', self.terms['0000036']['CHEBI'])


if __name__ == "__main__":
    unittest.main()
