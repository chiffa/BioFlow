'''
Created on Jun 15, 2013

@author: andrei
'''
import logging
from bulbs.neo4jserver import Graph as Neo4jGraph
import Neo4j_typeDec_new as DDT

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='dynamics_full.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

class Graph(Neo4jGraph):
    
    def __init__(self, config=None):
        super(Graph, self).__init__(config)
        
        #Annotations
        self.Location=self.build_proxy(DDT.Location)
        self.AnnotNode=self.build_proxy(DDT.AnnotNode)
        self.Originating_Organism=self.build_proxy(DDT.Originating_Organism)
        
        #Simple Compounds
        self.DNA=self.build_proxy(DDT.DNA)
        self.RNA=self.build_proxy(DDT.RNA)
        self.Protein=self.build_proxy(DDT.Protein)
        self.SmallMolecule=self.build_proxy(DDT.SmallMolecule)
        self.PhysicalEntity=self.build_proxy(DDT.PhysicalEntity)
        
        #And composite ones
        self.Complex=self.build_proxy(DDT.Complex)
        #That can be instantiated
        self.Instance=self.build_proxy(DDT.Instance)
        #By possible instantiating events
        self.ModificationFeature=self.build_proxy(DDT.ModificationFeature)
        
        #And that can be grouped as collections
        self.DNA_Collection=self.build_proxy(DDT.DNA_Collection)
        self.RNA_Collection=self.build_proxy(DDT.RNA_Collection)
        self.Protein_Collection=self.build_proxy(DDT.Protein_Collection)
        self.SmallMolecule_Collection=self.build_proxy(DDT.SmallMolecule_Collection)
        self.PhysicalEntity_Collection=self.build_proxy(DDT.PhysicalEntity_Collection)
        self.Complex_Collection=self.build_proxy(DDT.Complex_Collection)
        
        #And especially participate into reaction sets
        self.TemplateReaction=self.build_proxy(DDT.TemplateReaction)
        self.Degradation=self.build_proxy(DDT.Degradation)
        self.BiochemicalReaction=self.build_proxy(DDT.BiochemicalReaction)
        
        #Pointers from the Simple Compounds to their annotations
        self.is_localized=self.build_proxy(DDT.is_localized)
        self.is_annotated=self.build_proxy(DDT.is_annotated)
        self.is_originating_in_organism=self.build_proxy(DDT.is_originating_in_organism)
        
        #And from Complex Compounds to the simple Compounds they are made of
        self.is_part_of_complex=self.build_proxy(DDT.is_part_of_complex)
        
        #That can be instantiated
        self.is_instantiating=self.build_proxy(DDT.is_instantiating)
        #With Instantiators
        self.is_an_instantiator=self.is_a_possible_instance(DDT.is_an_instantiator)
        #Or belong to collections        
        self.is_part_of_collection=self.build_proxy(DDT.is_part_of_collection)
        
        #And contribute to reactions
        self.is_catalysant=self.build_proxy(DDT.is_catalysant)
        self.is_regulant=self.build_proxy(DDT.is_regulant)
        self.is_reaction_participant=self.build_proxy(DDT.is_reaction_participant)
        
        
#now, let's connect the graph and fill it with data from the etree parsing

DatabaseGraph=Graph()

# The idea is to do as little operations on the neo4j database, so we will be first loading everything in a system of dictionnaries and then
# flushing it alltogether to the neo4j database.

