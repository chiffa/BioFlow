"""
Created on Jul 2, 2013

@author: andrei
"""
from bulbs.neo4jserver import Graph as Neo4jGraph, Config
import neo4j_typeDec as DDL
from BioFlow.configs2 import neo4j_server
from BioFlow.Utils.LogManager import logger
import os

if neo4j_server != 'http://localhost:7474':
    neo4j_server_local = Config(neo4j_server + '/db/data/')
else:
    neo4j_server_local = None


class Graph(Neo4jGraph):
    """
    The interface with the neo4j graph database
    """

    def __init__(self, config=neo4j_server_local):
        super(Graph, self).__init__(config)
        
        # Annotations
        self.Location = self.build_proxy(DDL.Location)
        self.AnnotNode = self.build_proxy(DDL.AnnotNode)
        self.Originating_Organism = self.build_proxy(DDL.Originating_Organism)
        self.Pathway = self.build_proxy(DDL.Pathway)
        self.PathwayStep = self.build_proxy(DDL.Pathway_Step)
        self.GOTerm = self.build_proxy(DDL.GOTerm)
        self.UNIPORT = self.build_proxy(DDL.UNIPROT)
        
        # Simple Compounds
        self.DNA = self.build_proxy(DDL.DNA)
        self.RNA = self.build_proxy(DDL.RNA)
        self.Protein = self.build_proxy(DDL.Protein)
        self.SmallMolecule = self.build_proxy(DDL.SmallMolecule)
        self.PhysicalEntity = self.build_proxy(DDL.PhysicalEntity)
        
        # And composite ones
        self.Complex = self.build_proxy(DDL.Complex)
        # That can be instantiated
        self.Instance = self.build_proxy(DDL.Instance)
        # By possible instantiating events
        self.ModificationFeature = self.build_proxy(DDL.ModificationFeature)
        
        # And that can be grouped as collections
        self.DNA_Collection = self.build_proxy(DDL.DNA_Collection)
        self.RNA_Collection = self.build_proxy(DDL.RNA_Collection)
        self.Protein_Collection = self.build_proxy(DDL.Protein_Collection)
        self.SmallMolecule_Collection = self.build_proxy(DDL.SmallMolecule_Collection)
        self.PhysicalEntity_Collection = self.build_proxy(DDL.PhysicalEntity_Collection)
        self.Complex_Collection = self.build_proxy(DDL.Complex_Collection)
        
        # And especially participate into reaction sets
        self.TemplateReaction = self.build_proxy(DDL.TemplateReaction)
        self.Degradation = self.build_proxy(DDL.Degradation)
        self.BiochemicalReaction = self.build_proxy(DDL.BiochemicalReaction)
        
        # Pointers from the Simple Compounds to their annotations
        self.is_localized = self.build_proxy(DDL.is_localized)
        self.is_annotated = self.build_proxy(DDL.is_annotated)
        self.is_originating_in_organism = self.build_proxy(DDL.is_originating_in_organism)
        self.is_part_of_pathway = self.build_proxy(DDL.is_part_of_pathway)
        self.is_next_in_pathway = self.build_proxy(DDL.is_next_in_pathway)
        
        # And from Complex Compounds to the simple Compounds they are made of
        self.is_part_of_complex = self.build_proxy(DDL.is_part_of_complex)
        
        # That can be instantiated
        self.is_modified_to = self.build_proxy(DDL.is_modified_to)
        # With Instantiators
        self.is_able_to_modify = self.build_proxy(DDL.is_able_to_modify)
        # Or belong to collections
        self.is_part_of_collection = self.build_proxy(DDL.is_part_of_collection)
        
        # And contribute to reactions
        self.is_catalysant = self.build_proxy(DDL.is_catalysant)
        self.is_regulant = self.build_proxy(DDL.is_regulant)
        # regulates not a reaction, but a compound activity
        self.is_reaction_participant = self.build_proxy(DDL.is_reaction_participant)
        
        # GOAnnotationTypes:
        self.is_go_annotation = self.build_proxy(DDL.is_go_annotation)
        self.is_a_go = self.build_proxy(DDL.is_a_go)
        self.is_part_of_go = self.build_proxy(DDL.is_part_of_go)
        self.is_same = self.build_proxy(DDL.is_same)
        
        # Interacts physically:
        self.is_interacting = self.build_proxy(DDL.is_interacting)
        self.is_weakly_interacting = self.build_proxy(DDL.is_weakly_interacting)


# Yes, I know what goes below here is ugly and shouldn't be in the production part of the code

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
on_unittest = os.environ.get('UNITTESTING', None) == 'True'

if on_rtd or on_unittest:
    logger.info('GraphDeclarator mocks DB connection instead of actually connecting to it')

    from mock import Mock as MagicMock

    class Mock(MagicMock):

        @classmethod
        def __getattr__(cls, name):
            return Mock()

        @classmethod
        def __getitem__(cls, name):
            return Mock()

    DatabaseGraph = Mock()

else:
    logger.info('GraphDeclarator is connecting to a real DB')
    DatabaseGraph = Graph()
