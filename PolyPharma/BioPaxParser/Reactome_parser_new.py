'''
Created on Jun 15, 2013

@author: andrei
'''

from bulbs.neo4jserver import Graph as Neo4jGraph

import xml.etree.ElementTree as ET
import logging
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

tree = ET.parse('/home/andrei/UCSD/Parsing_Reactome/Homo sapiens.owl')
root = tree.getroot()

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

def zipDicts(dict1,dict2):
    '''
    performs a simple dictionarry assembly, adding up elements contained in 
    lists and throwing error in case non-list entries are zipped
    '''
    for key in dict2.keys():
        if key not in dict1.keys():
            dict1[key]=dict2[key]
        else:
            assert isinstance(dict2[key],(list,tuple))
            dict1[key]=dict1[key]+dict2[key]
    return dict1


################################################################3
# Simplest, pre-compression parses

BioSources={} #{ID:name}
BioSourcesXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}BioSource')
for single_bio_source in BioSourcesXml:
    key=single_bio_source.attrib.values()[0]
    for bio_source_ref_property in single_bio_source:
        if '}name' in bio_source_ref_property:
            BioSources[key]=bio_source_ref_property.text
            break

CellularLocations={} #{Id:name}
CellularLocationsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}CellularLocationVocabulary')
for single_cellular_location in CellularLocationsXml:
    key=single_cellular_location.attrib.values()[0]
    for loc_ref in single_cellular_location:
        if '}term' in loc_ref.tag:
            CellularLocations[key]=loc_ref.text


##################################################################
# a little bit more complicated parses

DnaRefs={} # {ID:{'name':[],'ENSEMBL':'','organism':''}}
DnaRefsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}ProteinReference')
for single_Dna_Ref in DnaRefsXml:
    key=single_Dna_Ref.attrib.values()[0]
    DnaRefs[key]={'name':[]}
    for dna_ref_property in single_Dna_Ref:
        if '}name' in dna_ref_property.tag:
            if 'ENSEMBL:' in dna_ref_property.text:
                line=dna_ref_property.text
                line=line.split()
                for word in line:
                    if 'ENSEMBL:' in word:
                        DnaRefs[key]['ENSEMBL']=word.split(':')[1]
            else:
                DnaRefs[key]['name'].append(dna_ref_property.text)
        if '}organism' in dna_ref_property.tag:
            DnaRefs[key]['organism']=BioSources[dna_ref_property.attrib.values()[0]]


####################################################################
# Core parses : parsing at the same time objects and their collections

def MetaParser_SecLoop(LocDic,local_property,CollectionMarker):
    if '}cellularLocation' in local_property.tag:
        LocDic['cellularLocation']=CellularLocations[local_property.attrib.values()[0]]
    if '}displayName' in local_property.tag:
        LocDic['displayName']=local_property.text
    if '}name' in local_property.tag:
        LocDic['references']['names'](local_property.text)
    if '}memberPhysicalEntity' in local_property.tag:
        CollectionMarker=True
        LocDic['collectionMemebers'].append(local_property.attrib.values()[0])


Dnas={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], references':{'names':[],...}}}
Dna_Collections={}
DnasXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Dna')
for single_Dna in DnasXml:
    key=single_Dna.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for dna_property in single_Dna:
        MetaParser_SecLoop(LocalDict, dna_property, Collection)
        if '}cellularLocation' in dna_property.tag:
            LocalDict['cellularLocation']=CellularLocations[dna_property.attrib.values()[0]]
        if '}displayName' in dna_property.tag:
            LocalDict['displayName']=dna_property.text
        if '}name' in dna_property.tag:
            LocalDict['references']['names'](dna_property.text)
        if '}entityReference' in dna_property.tag:
            LocalDict['references']=zipDicts(DnaRefs[dna_property.attrib.values()[0]], LocalDict['references'])
        if '}memberPhysicalEntity' in dna_property.tag:
            Collection=True
            LocalDict['collectionMemebers'].append(dna_property.attrib.values()[0])
    if Collection:
        Dna_Collections[key]=LocalDict
    else:
        Dnas[key]=LocalDict

PhysicalEntities={}
PhysicalEntity_Collections={}
PhysicalEntities_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity')
for single_PhysicalEntity in PhysicalEntities_Xml:
    key=single_PhysicalEntity.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for PhysicalEntity_property in single_PhysicalEntity:
        if '}cellularLocation' in PhysicalEntity_property.tag:
            LocalDict['cellularLocation']=CellularLocations[PhysicalEntity_property.attrib.values()[0]]
        if '}displayName' in PhysicalEntity_property.tag:
            LocalDict['displayName']=PhysicalEntity_property.text
        if '}name' in PhysicalEntity_property.tag:
            LocalDict['references']['names'](PhysicalEntity_property.text)
        if '}memberPhysicalEntity' in PhysicalEntity_property.tag:
            Collection=True
            LocalDict['collectionMemebers'].append(PhysicalEntity_property.attrib.values()[0])
    if Collection:
        PhysicalEntity_Collections[key]=LocalDict
    else:
        PhysicalEntities[key]=LocalDict





Complexes={} #{Id:{'cellularLocation':'', displayName:'', 'parts':[], collectionMembers':[], references':{'names':[],...}}}
Complex_Collections={}
ComplexXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Complex')
for single_Complex in ComplexXml:
    key=single_Complex.attrib.values()[0]
    LocalDict={'collectionMembers':[],'parts':[],'references':{'names':[]}}
    Collection=False
    for complex_property in single_Complex:
        if '}cellularLocation' in complex_property.tag:
            LocalDict['cellularLocation']=CellularLocations[complex_property .attrib.values()[0]]
        if '}displayName' in complex_property.tag:
            LocalDict['displayName']=complex_property.text
        if '}name' in complex_property .tag:
            LocalDict['references']['names'].append(complex_property.text)
        if '}memberPhysicalEntity' in complex_property.tag:
            Collection=True
            LocalDict['collectionMemebers'].append(complex_property.attrib.values()[0])
        if '}component' in complex_property.tag:
            LocalDict['parts'].append(complex_property.attrib.values()[0])
    if Collection:
        Complex_Collections[key]=LocalDict
    else:
        Complexes[key]=LocalDict



Catalysises={}#{ID:{Controller:'', Controlled:'', controlType:''}}
# TODO: first parse all thge controllers and keep adding the Flattened list of 
# catalyses they are participating to
# TODO: parse Modulations in a pre-compressed fashion and 
CatalysisesXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Catalysis')
for single_Catalysis in CatalysisesXml:
    key=single_Catalysis.atrrib.values()[0]
    for catalysis_property in single_Catalysis:
        if '}controlled' in catalysis_property.tag:
            Catalysises['controlled']=catalysis_property.attrib.values()[0]
        if '}controller' in catalysis_property.tag:
            Catalysises['controller']=catalysis_property.attrib.values()[0]
        if '}controlType' in catalysis_property.tag:
            Catalysises['ControlType']=catalysis_property.text

Modulations={} #{ID:{modulator, modulated} # This is essentially a compressed regulation of activity
# of the catalysts
ModulationsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Modulations')
for single_Modulation in ModulationsXml:
    key=single_Modulation.attrib.value()[0]
    for modulation_property in single_Modulation:
        if '}displayName' in modulation_property.tag:
            Modulations[key]['displayName']=modulation_property.text
        if '}controller' in modulation_property.tag:
            Modulations[key]['controller']=modulation_property.attrib.values()[0]
        if '}controlled' in modulation_property.tag:
            Modulations[key]['controlled']=Catalysises[modulation_property.attrib.values()[0]]['controller']
        if '}controlType' in modulation_property.tag:
            Modulations[key]['controlType']=modulation_property.text

