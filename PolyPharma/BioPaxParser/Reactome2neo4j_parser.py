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

UnificationXref={}
UnificationXrefs=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}UnificationXref')
for singlexref in UnificationXrefs:
    key=singlexref.attrib.values()[0]
    val=(singlexref[0].text, singlexref[1].text)
    UnificationXref[key]=val    
# Xref## : (dbName, ID within the DB)
# Ok, so for the proteins the Unification Xrefs are useless: they will have to be parsed from
# the list of human proteins before we can do any annotation on them. However the names are pretty 
# well defined

ProteinXref={}
ProteinXrefs=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}UnificationXref')
for single_prot_ref in ProteinXrefs:
    key=single_prot_ref.attrib.values()[0]
    val=[]
    for prot_ref_property in single_prot_ref:
        if '}name' in prot_ref_property.tag:
            val.append(prot_ref_property.text)
    ProteinXref[key]=val


DispName2Obj={}
DbNames=set()
Protein_Internal_Ref2Obj={}
Proteins=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Protein')
for single_prot in Proteins:
    LoadingObject={'xrefs':[], 'names':[], 'displayNames':[],'memberOf':[],'ProtRef':[]}
    for prot_property in single_prot:
        if '}xref' in prot_property.tag:
            LoadingObject['xrefs'].append(UnificationXref[prot_property.attrib.values()[0][1:]])
        if '}name' in prot_property.tag:
            LoadingObject['names'].append(prot_property.text)
        if '}displayName' in prot_property.tag:
            LoadingObject['displayNames'].append(prot_property.text)
#         if '}memberPhysicalEntity' in prot_property.tag:
#             LoadingObject['memberOf'].append(prot_property.attrib.values()[0])
        if '}entityReference' in prot_property.tag:
            LoadingObject['ProtRef'].append(prot_property.attrib.values()[0])
    
    for displayName in LoadingObject['displayNames']:
        if displayName not in DispName2Obj.keys():
            DispName2Obj[displayName]=[]
        DispName2Obj[displayName].append((LoadingObject['xrefs'], LoadingObject['names'],LoadingObject['memberOf'],LoadingObject['ProtRef']))
    Protein_Internal_Ref2Obj[single_prot.attrib.values()[0]]=LoadingObject

#Consider only the proteins that have a non-null external protein reference

for key in Protein_Internal_Ref2Obj.keys()[:10]:
    if Protein_Internal_Ref2Obj[key]['ProtRef']!=[]:
        print ProteinXref[Protein_Internal_Ref2Obj[key]['ProtRef'][0]]
#        No issues found: there are at most one prot.reference per protein
#         if len(Protein_Internal_Ref2Obj[key]['ProtRef'])>1:
#             print Protein_Internal_Ref2Obj[key]['ProtRef']

duplicates={}
for key in DispName2Obj.keys():
    if len(DispName2Obj[key])>1:
        if len(DispName2Obj[key]) not in duplicates.keys():
            duplicates[len(DispName2Obj[key])]=0
        duplicates[len(DispName2Obj[key])]+=1

print duplicates

a_n_plicate={}
for key in DispName2Obj.keys():
    if len(DispName2Obj[key])==8:
        a_n_plicate[key]=DispName2Obj[key]

for key in a_n_plicate.keys():
    print key 
    for lline in a_n_plicate[key]:
        print '\t', lline
    

# Now let's see if there are any conflicting names within the database            
# Up, there are quite a lot of them. Reasons:
#     - cellular location defines a protein, not it's instance
#     - Modification features (phosphorilation et Co)
#     - Fragments of chains within the protein structure (refered as fragment/domains for me)
# For me right now these are instances, under meta. Should I incorporate them?
# What we will do - we will invoke them all at once

# Whitelist the proteins that are adressed by the reaction.
# Problem: this will also change the instantiation methodology

# However the references to the Uniprot accession numbers are stored in the Protein Reference Field


Small_Molecules=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SmallMolecule')
Complex=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Complex')



SequenceMods=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
Cellular_Location=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}CellularLocationVocabulary')
Regulation=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Control')

#! Attention not to include the GO terms
Xref2=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}xref')
Xref1=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}RelationshipXref')
PathwaySteps=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PathwayStep')
PathwaySteps=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Pathway')

#Complex parsings:
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}FragmentFeature')
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReaction')
FragmentFeatures=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
PhysicalEntity=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity')
SeqInterval=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceInterval')
SeqSite=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceSite')
