import os
import unittest
from pprint import pprint
from BioFlow.utils.IO_Routines import dump_object, undump_object
from BioFlow.bio_db_parsers.gene_ontology_parser import GOTermsParser
from BioFlow.bio_db_parsers.new_unioprot_paser import UniProtParser


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


class UniprotParserTester(unittest.TestCase):

    up_to_parse = os.path.join(os.path.dirname(__file__), 'UT_examples/test_uniprot.dat')
    ref_parses = os.path.join(os.path.dirname(__file__), 'UT_examples/ref_up_parse.dmp')

    @classmethod
    def setUpClass(cls):
        parser_object = UniProtParser(['199310', '405955'])
        cls.uniprot_dict = parser_object.parse_uniprot(cls.up_to_parse)
        cls.acces_dict = parser_object.get_access_dicts()
        cls.ref_uniprot_dict, cls.ref_acces_dict = undump_object(cls.ref_parses)

    def test_total(self):  # TODO: in future, expand in smaller set of tests
        self.assertDictEqual(self.uniprot_dict, self.ref_uniprot_dict)
        self.assertDictEqual(self.acces_dict, self.ref_acces_dict)

if __name__ == "__main__":
    unittest.main()
