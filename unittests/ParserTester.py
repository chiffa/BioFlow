import os
import unittest
from pprint import pprint
from bioflow.utils.io_routines import dump_object, undump_object
from bioflow.bio_db_parsers.geneOntologyParser import GOTermsParser
from bioflow.bio_db_parsers.uniprotParser import UniProtParser
from bioflow.bio_db_parsers.reactomeParser import ReactomeParser


class GoParserTester(unittest.TestCase):

    ref_obo = os.path.join(os.path.dirname(__file__), 'UT_examples/test_go.obo')

    @classmethod
    def setUpClass(cls):
        cls.terms, cls.term_rels = GOTermsParser().parse_go_terms(cls.ref_obo)

    def test_proper_parsing(self):
        self.assertIn('0000001', list(self.terms.keys()))
        self.assertIn('0000002', list(self.terms.keys()))

    def test_namespaces(self):
        self.assertEqual('biological_process', self.terms['0000001']['namespace'])

    def test_names(self):
        self.assertEqual('mitochondrion inheritance', self.terms['0000001']['name'])

    def test_relations(self):
        self.assertIn(('0000001', 'is_a', '0048311'), self.term_rels)
        self.assertIn(('0000001', 'is_a', '0048308'), self.term_rels)
        self.assertIn(('0000002', 'is_a', '0007005'), self.term_rels)

    def test_obsolescence(self):
        self.assertNotIn('0000039', list(self.terms.keys()))

    def test_chebi_parsing(self):
        self.assertIn('CHEBI', list(self.terms['0000036'].keys()))
        self.assertIn('22221', self.terms['0000036']['CHEBI'])


class UniprotParserTester(unittest.TestCase):

    up_to_parse = os.path.join(os.path.dirname(__file__), 'UT_examples/test_uniprot.dat')
    ref_parses = os.path.join(os.path.dirname(__file__), 'UT_examples/ref_up_parse.dmp')

    @classmethod
    def setUpClass(cls):
        parser_object = UniProtParser(['199310', '405955'])
        cls.uniprot_dict = parser_object.parse_uniprot(cls.up_to_parse)
        cls.acces_dict = parser_object.get_access_dicts()
        # dump_object(cls.ref_parses, (cls.uniprot_dict, cls.acces_dict))
        cls.ref_uniprot_dict, cls.ref_acces_dict = undump_object(cls.ref_parses)

    def test_total(self):
        self.assertDictEqual(self.uniprot_dict, self.ref_uniprot_dict)
        self.assertDictEqual(self.acces_dict, self.ref_acces_dict)


# class ReactomeParseTester(unittest.TestCase):
#     # Developers, please notice that for the real modification you would need to load the real
#     # reactome and run tests against it's initial dump first.
#     reactome_to_parse = os.path.join(os.path.dirname(__file__), 'UT_examples/reactome_extract.owl')
#     ref_parse = os.path.join(os.path.dirname(__file__), 'UT_examples/ref_reactome_parse.dmp')
#
#     @classmethod
#     def setUpClass(cls):
#         cls.actual_parser = ReactomeParser(cls.reactome_to_parse)
#         cls.actual_parser.parse_all()
#         # print(cls.actual_parser.parsed)
#         # dump_object(cls.ref_parse, cls.actual_parser)
#         cls.ref_parser = undump_object(cls.ref_parse)
#         # to parse
#         cls.maxDiff = None
#
#     def test_BioSources(self):
#         self.assertDictEqual(self.actual_parser.BioSources,
#                              self.ref_parser.BioSources)
#
#     def test_CellularLocations(self):
#         self.assertDictEqual(self.actual_parser.CellularLocations,
#                              self.ref_parser.CellularLocations)
#
#     def test_SeqModVoc(self):
#         self.assertDictEqual(self.actual_parser.SeqModVoc,
#                              self.ref_parser.SeqModVoc)
#
#     def test_SeqSite(self):
#         self.assertDictEqual(self.actual_parser.SeqSite,
#                              self.ref_parser.SeqSite)
#
#     def test_DnaRefs(self):
#         self.assertDictEqual(self.actual_parser.DnaRefs,
#                              self.ref_parser.DnaRefs)
#
#     def test_RnaRefs(self):
#         self.assertCountEqual(self.actual_parser.RnaRefs,
#                              self.ref_parser.RnaRefs)
#
#     def test_SmallMoleculeRefs(self):
#         self.assertCountEqual(self.actual_parser.SmallMoleculeRefs,
#                              self.ref_parser.SmallMoleculeRefs)
#
#     def test_ProteinRefs(self):
#         self.assertCountEqual(self.actual_parser.ProteinRefs,
#                              self.ref_parser.ProteinRefs)
#
#     def test_ModificationFeatures(self):
#         self.assertDictEqual(self.actual_parser.ModificationFeatures,
#                              self.ref_parser.ModificationFeatures)
#
#     def test_Dna(self):
#         self.assertDictEqual(self.actual_parser.Dnas,
#                              self.ref_parser.Dnas)
#         self.assertDictEqual(self.actual_parser.Dna_Collections,
#                              self.ref_parser.Dna_Collections)
#
#     def test_Rna(self):
#         self.assertCountEqual(self.actual_parser.Rnas,
#                              self.ref_parser.Rnas)
#         self.assertCountEqual(self.actual_parser.Rna_Collections,
#                              self.ref_parser.Rna_Collections)
#
#     def test_SmallMolecules(self):
#         self.assertCountEqual(self.actual_parser.SmallMolecules,
#                              self.ref_parser.SmallMolecules)
#         self.assertCountEqual(self.actual_parser.SmallMolecule_Collections,
#                              self.ref_parser.SmallMolecule_Collections)
#
#     def test_Proteins(self):
#         self.assertCountEqual(self.actual_parser.Proteins,
#                              self.ref_parser.Proteins)
#         self.assertCountEqual(self.actual_parser.Protein_Collections,
#                              self.ref_parser.Protein_Collections)
#
#     def test_PhysicalEntities(self):
#         self.assertDictEqual(self.actual_parser.PhysicalEntities,
#                              self.ref_parser.PhysicalEntities)
#         self.assertDictEqual(self.actual_parser.PhysicalEntity_Collections,
#                              self.ref_parser.PhysicalEntity_Collections)
#
#     def test_Complexes(self):
#         self.assertDictEqual(self.actual_parser.Complexes,
#                              self.ref_parser.Complexes)
#         self.assertDictEqual(self.actual_parser.Complex_Collections,
#                              self.ref_parser.Complex_Collections)
#
#     def test_TemplateReactions(self):
#         self.assertDictEqual(self.actual_parser.TemplateReactions,
#                              self.ref_parser.TemplateReactions)
#
#     def test_Degradations(self):
#         self.assertDictEqual(self.actual_parser.Degradations,
#                              self.ref_parser.Degradations)
#
#     def test_BiochemicalReactions(self):
#         self.assertDictEqual(self.actual_parser.BiochemicalReactions,
#                              self.ref_parser.BiochemicalReactions)
#
#     def test_Catalysises(self):
#         self.assertDictEqual(self.actual_parser.Catalysises,
#                              self.ref_parser.Catalysises)
#
#     def test_Modulation(self):
#         self.assertDictEqual(self.actual_parser.Modulations,
#                              self.ref_parser.Modulations)
#
#     def test_Pathways(self):
#         self.assertDictEqual(self.actual_parser.Pathways,
#                              self.ref_parser.Pathways)
#
#     def test_PathwaySteps(self):
#         self.assertDictEqual(self.actual_parser.PathwaySteps,
#                              self.ref_parser.PathwaySteps)


if __name__ == "__main__":
    unittest.main()
