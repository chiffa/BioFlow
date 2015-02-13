'''
Created on Jul 2, 2013

@author: andrei
'''

if __name__ == "__main__" and __package__ is None:
    __package__ = "PolyPharma.neo4j_Declarations"

from bulbs.neo4jserver import Graph as Neo4jGraph, Config
import neo4j_typeDec as DDT
from PolyPharma.configs import neo4j_server

if neo4j_server != 'http://localhost:7474':
    neo4j_server_local = Config(neo4j_server+ '/db/data/')
else: neo4j_server_local = None

# noinspection PyTypeChecker
class Graph(Neo4jGraph):
    '''
    The interface with the neo4j graph database
    '''

    def __init__(self, config=neo4j_server_local):
        super(Graph, self).__init__(config)
        
        #Annotations
        self.Location = self.build_proxy(DDT.Location)
        self.AnnotNode = self.build_proxy(DDT.AnnotNode)
        self.Originating_Organism = self.build_proxy(DDT.Originating_Organism)
        self.Pathway = self.build_proxy(DDT.Pathway)
        self.PathwayStep = self.build_proxy(DDT.Pathway_Step)
        self.GOTerm = self.build_proxy(DDT.GOTerm)
        self.UNIPORT = self.build_proxy(DDT.UNIPROT)
        
        #Simple Compounds
        self.DNA = self.build_proxy(DDT.DNA)
        self.RNA = self.build_proxy(DDT.RNA)
        self.Protein = self.build_proxy(DDT.Protein)
        self.SmallMolecule = self.build_proxy(DDT.SmallMolecule)
        self.PhysicalEntity = self.build_proxy(DDT.PhysicalEntity)
        
        #And composite ones
        self.Complex = self.build_proxy(DDT.Complex)
        #That can be instantiated
        self.Instance = self.build_proxy(DDT.Instance)
        #By possible instantiating events
        self.ModificationFeature = self.build_proxy(DDT.ModificationFeature)
        
        #And that can be grouped as collections
        self.DNA_Collection = self.build_proxy(DDT.DNA_Collection)
        self.RNA_Collection = self.build_proxy(DDT.RNA_Collection)
        self.Protein_Collection = self.build_proxy(DDT.Protein_Collection)
        self.SmallMolecule_Collection = self.build_proxy(DDT.SmallMolecule_Collection)
        self.PhysicalEntity_Collection = self.build_proxy(DDT.PhysicalEntity_Collection)
        self.Complex_Collection = self.build_proxy(DDT.Complex_Collection)
        
        #And especially participate into reaction sets
        self.TemplateReaction = self.build_proxy(DDT.TemplateReaction)
        self.Degradation = self.build_proxy(DDT.Degradation)
        self.BiochemicalReaction = self.build_proxy(DDT.BiochemicalReaction)
        
        #Pointers from the Simple Compounds to their annotations
        self.is_localized = self.build_proxy(DDT.is_localized)
        self.is_annotated = self.build_proxy(DDT.is_annotated)
        self.is_originating_in_organism = self.build_proxy(DDT.is_originating_in_organism)
        self.is_part_of_pathway = self.build_proxy(DDT.is_part_of_pathway)
        self.is_next_in_pathway = self.build_proxy(DDT.is_next_in_pathway)
        
        #And from Complex Compounds to the simple Compounds they are made of
        self.is_part_of_complex = self.build_proxy(DDT.is_part_of_complex)
        
        #That can be instantiated
        self.is_modified_to = self.build_proxy(DDT.is_modified_to)
        #With Instantiators
        self.is_able_to_modify = self.build_proxy(DDT.is_able_to_modify)
        #Or belong to collections        
        self.is_part_of_collection = self.build_proxy(DDT.is_part_of_collection)
        
        #And contribute to reactions
        self.is_catalysant = self.build_proxy(DDT.is_catalysant)
        self.is_regulant = self.build_proxy(DDT.is_regulant)  # regulates not a reaction, but a compound activity
        self.is_reaction_participant = self.build_proxy(DDT.is_reaction_participant)
        
        #GOAnnotationTypes:
        self.is_go_annotation = self.build_proxy(DDT.is_go_annotation)
        self.is_a_go = self.build_proxy(DDT.is_a_go)
        self.is_part_of_go = self.build_proxy(DDT.is_part_of_go)
        self.is_same = self.build_proxy(DDT.is_same)
        
        # Interacts physically:
        self.is_interacting = self.build_proxy(DDT.is_interacting)
        self.is_weakly_interacting = self.build_proxy(DDT.is_weakly_interacting)
        

DatabaseGraph = Graph()