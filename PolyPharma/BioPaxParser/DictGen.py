'''
Created on Jun 17, 2013

@author: andrei
#Put it into the a Reactome parser package later on
'''
# TODO: we might want to parse the traceability of the all the compouunds and link by adding
# the xref parsed information to them

#TODO: delete the '#' for the references added from within the properties

import logging
import xml.etree.ElementTree as ET

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

SeqModVoc={} #{Id:name}
SeqModVocXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
for single_SeqMod in SeqModVocXml:
    key=single_SeqMod.attrib.values()[0]
    for mod_ref in single_SeqMod:
        if '}term' in mod_ref.tag:
            SeqModVoc[key]=mod_ref.text

SeqSite={} #{Id:Position}
SeqSiteXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
for single_SeqSite in SeqSiteXml:
    key=single_SeqSite.attrib.values()[0]
    for site_ref in single_SeqSite:
        if '}sequencePosition' in site_ref.tag:
            SeqSite[key]=site_ref.text

##################################################################
# a little bit more complicated reference parses

DnaRefs={} # {ID:{'name':[],'ENSEMBL':'','organism':''}}
DnaRefsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}DnaReference')
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


RnaRefs={} # {ID:{'name':[],'ENSEMBL':'', 'miRBase':'', EMBL:'', s'organism':''}}
RnaRefsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}RnaReference')
for single_Rna_Ref in RnaRefsXml:
    key=single_Rna_Ref.attrib.values()[0]
    RnaRefs[key]={'name':[]}
    for rna_ref_property in single_Rna_Ref:
        if '}name' in rna_ref_property.tag:
            if 'ENSEMBL:' in rna_ref_property.text:
                line=rna_ref_property.text
                line=line.split()
                for word in line:
                    if 'ENSEMBL:' in word:
                        RnaRefs[key]['ENSEMBL']=word.split(':')[1]
            if 'miRBase:' in rna_ref_property.text:
                line=rna_ref_property.text
                line=line.split()
                for word in line:
                    if 'miRBase:' in word:
                        RnaRefs[key]['miRBase']=word.split(':')[1]
            if 'EMBL:' in rna_ref_property.text:
                line=rna_ref_property.text
                line=line.split()
                for word in line:
                    if 'EMBL:' in word:
                        RnaRefs[key]['EMBL']=word.split(':')[1]
            else:
                RnaRefs[key]['name'].append(rna_ref_property.text)
        if '}organism' in rna_ref_property.tag:
            RnaRefs[key]['organism']=BioSources[rna_ref_property.attrib.values()[0]]

SmallMoleculeRefs={} # ID:{'name':[],'ChEBI':'','organism':''}}
SmallMoleculeRefsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SmallMoleculeReference')
for single_SmallMolecule_Ref in SmallMoleculeRefsXml:
    key=single_SmallMolecule_Ref.attrib.values()[0]
    SmallMoleculeRefs[key]={'name':[]}
    for SmallMolecule_ref_property in single_SmallMolecule_Ref:
        if '}name' in SmallMolecule_ref_property.tag:
            if 'ChEBI:' in SmallMolecule_ref_property.text:
                line=SmallMolecule_ref_property.text
                line=line.split()
                for word in line:
                    if 'ChEBI:' in word:
                        SmallMoleculeRefs[key]['ChEBI']=word.split(':')[1].strip(']')
            else:
                SmallMoleculeRefs[key]['name'].append(SmallMolecule_ref_property.text)
        if '}organism' in SmallMolecule_ref_property.tag:
            SmallMoleculeRefs[key]['organism']=BioSources[SmallMolecule_ref_property.attrib.values()[0]]



ProteinRefs={} # {ID:{'name':[],'ENSEMBL':'','organism':''}}
ProteinRefsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}ProteinReference')
for single_Protein_Ref in ProteinRefsXml:
    key=single_Protein_Ref.attrib.values()[0]
    ProteinRefs[key]={'name':[]}
    for protein_ref_property in single_Protein_Ref:
        if '}name' in protein_ref_property.tag:
            if 'UniProt:' in protein_ref_property.text:
                line=protein_ref_property.text
                line=line.split()
                for word in line:
                    if 'UniProt:' in word:
                        ProteinRefs[key]['UniProt']=word.split(':')[1]
            else:
                ProteinRefs[key]['name'].append(protein_ref_property.text)
        if '}organism' in protein_ref_property.tag:
            ProteinRefs[key]['organism']=BioSources[protein_ref_property.attrib.values()[0]]


ModificationFeatures={} # {ID:{'location':'','modification':''}}
ModificationFeatureXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}ModificationFeature')
for single_ModificationFeature in ModificationFeatureXml:
    key=single_ModificationFeature.attrib.values()[0]
    ModificationFeatures[key]={}
    for modification_property in single_ModificationFeature:
        if '}featureLocation' in modification_property.tag:
            ModificationFeatures[key]['location']=SeqSite[modification_property.attrib.values()[0]]
        if '}modificationType' in modification_property.tag:
            ModificationFeatures[key]['modification']=SeqModVoc[modification_property.attrib.values()[0]]
            

####################################################################
# Core parses : parsing at the same time objects and their collections

def MetaParser_SecLoop(LocDic,local_property,CollectionMarker):
    '''
    this part is the same for all the Meta objects, so let's save some space
    '''
    
    if '}cellularLocation' in local_property.tag:
        LocDic['cellularLocation']=CellularLocations[local_property.attrib.values()[0]]
    if '}displayName' in local_property.tag:
        LocDic['displayName']=local_property.text
    if '}name' in local_property.tag:
        LocDic['references']['names'](local_property.text)
    if '}memberPhysicalEntity' in local_property.tag:
        CollectionMarker=True
        LocDic['collectionMemebers'].append(local_property.attrib.values()[0])
    return CollectionMarker

Dnas={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], references':{'names':[],...}}}
Dna_Collections={}
DnasXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Dna')
for single_Dna in DnasXml:
    key=single_Dna.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for dna_property in single_Dna:
        Collection=MetaParser_SecLoop(LocalDict, dna_property, Collection)
        if '}entityReference' in dna_property.tag:
            LocalDict['references']=zipDicts(DnaRefs[dna_property.attrib.values()[0]], LocalDict['references'])
    if Collection:
        Dna_Collections[key]=LocalDict
    else:
        Dnas[key]=LocalDict

Rnas={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], references':{'names':[],...}}}
Rna_Collections={}
RnasXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Rna')
for single_Rna in RnasXml:
    key=single_Rna.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for Rna_property in single_Rna:
        Collection=MetaParser_SecLoop(LocalDict, Rna_property, Collection)
        if '}entityReference' in Rna_property.tag:
            LocalDict['references']=zipDicts(RnaRefs[Rna_property.attrib.values()[0]], LocalDict['references'])
        Rna_Collections[key]=LocalDict
    else:
        Rnas[key]=LocalDict

SmallMolecules={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], references':{'names':[],...}}}
SmallMolecule_Collections={}
SmallMolecules_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SmallMolecules')
for single_SmallMolecule in SmallMolecules_Xml:
    key=single_SmallMolecule.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for SmallMolecule_property in single_SmallMolecule:
        Collection=MetaParser_SecLoop(LocalDict, SmallMolecule_property, Collection)
        if '}entityReference' in Rna_property.tag:
            LocalDict['references']=zipDicts(SmallMoleculeRefs[SmallMolecule_property.attrib.values()[0]], LocalDict['references'])
    if Collection:
        SmallMolecule_Collections[key]=LocalDict
    else:
        SmallMolecules[key]=LocalDict

PhysicalEntities={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], references':{'names':[],...}}}
PhysicalEntity_Collections={}
PhysicalEntities_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity')
for single_PhysicalEntity in PhysicalEntities_Xml:
    key=single_PhysicalEntity.attrib.values()[0]
    LocalDict={'collectionMembers':[],'references':{'names':[]}}
    Collection=False
    for PhysicalEntity_property in single_PhysicalEntity:
        Collection=MetaParser_SecLoop(LocalDict, PhysicalEntity_property, Collection)
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
        Collection=MetaParser_SecLoop(LocalDict, PhysicalEntity_property, Collection)
        if '}component' in complex_property.tag:
            LocalDict['parts'].append(complex_property.attrib.values()[0])
    if Collection:
        Complex_Collections[key]=LocalDict
    else:
        Complexes[key]=LocalDict

TemplateReactions={}
TemplateReactionXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
for single_TemplateReaction in TemplateReactionXml:
    key=single_Complex.attrib.values()[0]
    LocalDict={'references':{'names':[]}}
    for TemplateReaction_property in single_Complex:
        if '}product' in complex_property.tag:
            LocalDict['product']=TemplateReaction_property.attrib.values()[0]
        if '}displayName' in TemplateReaction_property.tag:
            LocalDict['displayName']=TemplateReaction_property.text
        if '}name' in TemplateReaction_property.tag:
            LocalDict['references']['names']=TemplateReaction_property.text
        TemplateReactions[key]=LocalDict

Degradation={}
DegradationXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
for single_Degradation in DegradationXml:
    key=single_Complex.attrib.values()[0]
    LocalDict={'references':{'names':[]}}
    for TemplateReaction_property in single_Complex:
        if '}left' in complex_property.tag:
            LocalDict['degraded']=TemplateReaction_property.attrib.values()[0]
        if '}displayName' in TemplateReaction_property.tag:
            LocalDict['displayName']=TemplateReaction_property.text
        if '}eCNumber' in TemplateReaction_property.tag:
            LocalDict['references']['eCNumber']=TemplateReaction_property.attrib.values()[0]
        TemplateReactions[key]=LocalDict

BiochemicalReactions={}
BiochemicalReactionXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
for single_BiochemicalReaction in BiochemicalReactionXml:
    key=single_BiochemicalReaction.attrib.values()[0]
    LocalDict={'right':[],'left':[],'references':{'names':[]}}
    for BiochemicalReaction_property in single_BiochemicalReaction:
        if '}left' in BiochemicalReaction_property.tag:
            LocalDict['left']=BiochemicalReaction_property.attrib.values()[0]
        if '}right' in BiochemicalReaction_property.tag:
            LocalDict['right']=BiochemicalReaction_property.attrib.values()[0]
        if '}displayName' in BiochemicalReaction_property.tag:
            LocalDict['displayName']=BiochemicalReaction_property.text
        if '}eCNumber' in BiochemicalReaction_property.tag:
            LocalDict['references']['eCNumber']=BiochemicalReaction_property.attrib.values()[0]
        TemplateReactions[key]=LocalDict

Catalysises={}#{ID:{Controller:'', Controlled:'', controlType:''}}
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
TemplateReactRegulationXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
for single_TRR in TemplateReactRegulationXml:
    key=single_Catalysis.atrrib.values()[0]
    for TRR_property in single_TRR:
        if '}controlled' in TRR_property.tag:
            Catalysises['controlled']=TRR_property.attrib.values()[0]
        if '}controller' in TRR_property.tag:
            Catalysises['controller']=TRR_property.attrib.values()[0]
        if '}controlType' in TRR_property.tag:
            Catalysises['ControlType']=TRR_property.text
ControlXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Control')
for single_Control in ControlXml:
    key=single_Control.attrib.values()[0]
    for control_property in single_Control:
        if '}displayName' in control_property.tag:
            Catalysises['displayName']=control_property.text
        if '}controlled' in control_property.tag:
            Catalysises['controlled']=control_property.attrib.values()[0]
        if '}controller' in control_property.tag:
            Catalysises['controller']=control_property.attrib.values()[0]
        if '}controlType' in control_property.tag and 'BiochemicalReaction' in control_property.attrib.values()[0]:
            Catalysises['ControlType']=control_property.text


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


