'''
Created on May 30, 2013

@author: andrei
'''

from bulbs.neo4jserver import Graph as Neo4jGraph

import xml.etree.ElementTree as ET
import logging
import Neo4j_Type_Declaration as DDT

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

tree = ET.parse('/home/andrei/UCSD/Parsing_Reactome/Homo sapiens.owl')
root = tree.getroot()


class Graph(Neo4jGraph):
    
    def __init__(self, config=None):
        super(Graph, self).__init__(config)
        
        # Node Proxies
        self.Protein = self.build_proxy(DDT.Meta_Protein)
        self.Chromosome = self.build_proxy(DDT.Meta_Chromosome)
        self.Rna = self.build_proxy(DDT.Meta_RNA)
        self.Complex=self.build_proxy(DDT.Meta_Complex)
        self.Small_molecule=self.build_proxy(DDT.Meta_SmallMolecule)
        
        self.Domain=self.build_proxy(DDT.Fragment_Domain)
        self.Dna=self.build_proxy(DDT.Fragment_DNA)
        self.Evosite=self.build_proxy(DDT.Fragment_EvoDefSite)
        self.Localization=self.build_proxy(DDT.Instantiation_Localization)
        
        self.Reaction=self.build_proxy(DDT.Reaction)
        
        self.Database=self.build_proxy(DDT.Annot_Database)
        self.Location=self.build_proxy(DDT.Annot_Location_type)
        self.Reaction_type=self.build_proxy(DDT.Annot_Reaction_type)
        self.GO_Term=self.build_proxy(DDT.Annot_GOTerm)
        self.Common_Name=self.build_proxy(DDT.Annot_CommonName)
        self.Pathway=self.build_proxy(DDT.Annot_Pathway)
        

        # Relationship Proxies
        self.rel = self.build_proxy(DDT.toA_CommonName)
        self.type = self.build_proxy(DDT.ToA_Typing)
        self.encodes = self.build_proxy(DDT.M2M_Encodes)
        self.maps_to = self.build_proxy(DDT.M2F_Maps_To)
        self.annotates = self.build_proxy(DDT.M_F2A_Annotates)
        self.xref = self.build_proxy(DDT.M_F_I2A_Xref)
        self.participates_in = self.build_proxy(DDT.M2R_Participates_in_reaction)
        self.regulates = self.build_proxy(DDT.M2R_Regulates)
        self.named = self.build_proxy(DDT.toA_CommonName)
        self.part_of = self.build_proxy(DDT.M_F2M_Part_of)
        self.in_pathway = self.build_proxy(DDT.R_A2A_InPatway)
        self.next_step = self.build_proxy(DDT.R2R_NextStepInPathway)
        
#now, let's connect the graph and fill it with data from the etree parsing

DatabaseGraph=Graph()

SequenceMods=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
Cellular_Location=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}CellularLocationVocabulary')
Small_Molecules=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SmallMolecule')
Complex=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Complex')
'{http://www.biopax.org/release/biopax-level3.owl#}Control'
'{http://www.biopax.org/release/biopax-level3.owl#}Protein'
#! Attention not to include the GO terms
'{http://www.biopax.org/release/biopax-level3.owl#}RelationshipXref'
'{http://www.biopax.org/release/biopax-level3.owl#}PathwayStep'
'{http://www.biopax.org/release/biopax-level3.owl#}Pathway'

#Complex parsings:
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}FragmentFeature')
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReaction')
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
'{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity'
'{http://www.biopax.org/release/biopax-level3.owl#}SequenceInterval'
'{http://www.biopax.org/release/biopax-level3.owl#}SequenceSite'
