'''
Created on Jun 17, 2013

@author: andrei
#Put it into the a Reactome parser package later on
'''

import logging
import xml.etree.ElementTree as ET
from random import shuffle
import PolyPharma.configs as conf

# TODO: refactor as a class

####################################################################################
#
# Logger behavior definition (most of the time it fails to function due to the)
# logs collision with neo4j
#
####################################################################################

# TODO: export logs location to the configs file

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='../logs/dynamics_full_2.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)


####################################################################################
#
# We define here the containers for all the values we are going to be interested in
#
####################################################################################

BioSources={} #{ID:name}
CellularLocations={} #{Id:name}
SeqModVoc={} #{Id:name}
SeqSite={} #{Id:Position}
Pathways={} #{ID:{'displayName':'','pathwayComponent':'','name':[], 'PathwayStep':[]}}
PathwaySteps={} #{ID:{'stepProcess':[],'nextStep':[]}}

DnaRefs={} # {ID:{'name':[],'ENSEMBL':'','organism':''}}
RnaRefs={} # {ID:{'name':[],'ENSEMBL':'', 'miRBase':'', EMBL:'', s'organism':''}}
SmallMoleculeRefs={} # ID:{'name':[],'ChEBI':'','organism':''}}
ProteinRefs={} # {ID:{'name':[],'UniProt':'','organism':''}}
ModificationFeatures={} # {ID:{'location':'','modification':''}}

Dnas={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{'names':[],...}}}
Dna_Collections={}
Rnas={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{'names':[],...}}}
Rna_Collections={}
SmallMolecules={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{'names':[],...}}}
SmallMolecule_Collections={}
Proteins={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{'names':[],...}}}
Protein_Collections={}
PhysicalEntities={} #{Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{'names':[],...}}}
PhysicalEntity_Collections={}
Complexes={} #{Id:{'cellularLocation':'', displayName:'', 'parts':[], collectionMembers':[], 'references':{'names':[],...}}}
Complex_Collections={}

TemplateReactions={} # {ID:{'product':'','displayName':'','references':{'names':[],...}}}
Degradations={}# {ID:{'product':'','displayName':'','references':{'eCNumber':[],...}}}
BiochemicalReactions={} # {ID:{'left':[],'right':[],'displayName':'','references':{'eCNumber':[],...}}}
Catalysises={}#{ID:{Controller:'', Controlled:'', controlType:''}}

Modulations={} #{ID:{modulator, modulated} # This is essentially a compressed regulation of activity of the catalysts


def shDic(dico, entries=10):
    """
    A supporting debug function.
    Extracts a random subset of a dictionary and prints it
    """
    EntryList = dico.items()
    shuffle(EntryList)
    print 'length:', len(EntryList)
    for elt in EntryList[0:entries]:
        print elt[0],'\t| ',elt[1]       

def zipDicts(dict1,dict2):
    '''
    performs a simple dictionarry assembly, adding up elements contained in 
    lists and throwing error in case non-list entries are zipped
    '''
    for key in dict2.keys():
        if key not in dict1.keys():
            dict1[key] = dict2[key]
        else:
            assert isinstance(dict2[key],(list,tuple))
            dict1[key]=dict1[key]+dict2[key]
    return dict1

####################################################################################
#
# Simplest, pre-compression parses
#
####################################################################################

def parse_BioSource(root):
    BioSourcesXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}BioSource')
    for single_bio_source in BioSourcesXml:
        key=single_bio_source.attrib.values()[0]
        for bio_source_ref_property in single_bio_source:
            if '}name' in bio_source_ref_property.tag:
                BioSources[key]=bio_source_ref_property.text
                break
    
def parse_CellularLocations(root):
    CellularLocationsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}CellularLocationVocabulary')
    for single_cellular_location in CellularLocationsXml:
        key=single_cellular_location.attrib.values()[0]
        for loc_ref in single_cellular_location:
            if '}term' in loc_ref.tag:
                CellularLocations[key]=loc_ref.text
    
def parse_SeqModVoc(root):
    SeqModVocXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
    for single_SeqMod in SeqModVocXml:
        key=single_SeqMod.attrib.values()[0]
        for mod_ref in single_SeqMod:
            if '}term' in mod_ref.tag:
                SeqModVoc[key]=mod_ref.text

def parse_SeqSite(root):
    SeqSiteXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SequenceSite')
    for single_SeqSite in SeqSiteXml:
        key=single_SeqSite.attrib.values()[0]
        for site_ref in single_SeqSite:
            if '}sequencePosition' in site_ref.tag:
                SeqSite[key]=site_ref.text

####################################################################################
#
# a little bit more complicated reference parses
#
####################################################################################

def parse_DnaRefs(root):
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
            if '}organism' in dna_ref_property.tag and dna_ref_property.attrib.values()[0][1:]!='BioSource1':
                DnaRefs[key]['organism']=BioSources[dna_ref_property.attrib.values()[0][1:]]
    
def parse_RnaRefs(root):
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
                else: 
                    if 'miRBase:' in rna_ref_property.text:
                        line=rna_ref_property.text
                        line=line.split()
                        for word in line:
                            if 'miRBase:' in word:
                                RnaRefs[key]['miRBase']=word.split(':')[1]
                    else: 
                        if 'EMBL:' in rna_ref_property.text:
                            line=rna_ref_property.text
                            line=line.split()
                            for word in line:
                                if 'EMBL:' in word:
                                    RnaRefs[key]['EMBL']=word.split(':')[1]
                        else:
                            RnaRefs[key]['name'].append(rna_ref_property.text)
            if '}organism' in rna_ref_property.tag and rna_ref_property.attrib.values()[0][1:]!='BioSource1':
                RnaRefs[key]['organism']=BioSources[rna_ref_property.attrib.values()[0][1:]]
    
def parse_SmallMoleculeRefs(root):
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
            if '}organism' in SmallMolecule_ref_property.tag and SmallMolecule_ref_property.attrib.values()[0][1:]!='BioSource1':
                SmallMoleculeRefs[key]['organism']=BioSources[SmallMolecule_ref_property.attrib.values()[0][1:]]
    
def parse_ProteinRefs(root):
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
            if '}organism' in protein_ref_property.tag and protein_ref_property.attrib.values()[0][1:]!='BioSource1':
                ProteinRefs[key]['organism']=BioSources[protein_ref_property.attrib.values()[0][1:]]
    
def parse_ModificationFeatures(root):
    ModificationFeatureXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}ModificationFeature')
    for single_ModificationFeature in ModificationFeatureXml:
        key=single_ModificationFeature.attrib.values()[0]
        ModificationFeatures[key]={'ID':key}
        for modification_property in single_ModificationFeature:
            if '}featureLocation' in modification_property.tag:
                ModificationFeatures[key]['location']=SeqSite[modification_property.attrib.values()[0][1:]]
            if '}modificationType' in modification_property.tag:
                ModificationFeatures[key]['modification']=SeqModVoc[modification_property.attrib.values()[0][1:]]

####################################################################################
#
# Core parses : parsing at the same time objects and their collections
#
####################################################################################


def MetaParser_SecLoop(LocDic,local_property,CollectionMarker):
    '''
    this part is the same for all the Meta objects, so let's save some space
    '''
    
    if '}cellularLocation' in local_property.tag:
        LocDic['cellularLocation']=local_property.attrib.values()[0][1:]
    if '}displayName' in local_property.tag:
        LocDic['displayName']=local_property.text
    if '}name' in local_property.tag:
        LocDic['references']['name'].append(local_property.text)
    if '}memberPhysicalEntity' in local_property.tag:
        CollectionMarker=True
        LocDic['collectionMembers'].append(local_property.attrib.values()[0][1:])
    return CollectionMarker

def parse_Dnas(root):
    Dnas_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Dna')
    for single_Dna in Dnas_Xml:
        key=single_Dna.attrib.values()[0]
        LocalDict={'collectionMembers':[],'references':{'name':[]}}
        Collection=False
        for dna_property in single_Dna:
            Collection=MetaParser_SecLoop(LocalDict, dna_property, Collection)
            if '}entityReference' in dna_property.tag:
                LocalDict['references']=zipDicts(DnaRefs[dna_property.attrib.values()[0][1:]], LocalDict['references'])
        if Collection:
            Dna_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            Dnas[key]=LocalDict
    
def parse_Rnas(root):
    Rnas_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Rna')
    for single_Rna in Rnas_Xml:
        key=single_Rna.attrib.values()[0]
        LocalDict={'collectionMembers':[],'references':{'name':[]}}
        Collection=False
        for Rna_property in single_Rna:
            Collection=MetaParser_SecLoop(LocalDict, Rna_property, Collection)
            if '}entityReference' in Rna_property.tag:
                LocalDict['references']=zipDicts(RnaRefs[Rna_property.attrib.values()[0][1:]], LocalDict['references'])
        if Collection:
            Rna_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            Rnas[key]=LocalDict
    
def parse_SmallMolecules(root):
    SmallMolecules_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}SmallMolecule')
    for single_SmallMolecule in SmallMolecules_Xml:
        key=single_SmallMolecule.attrib.values()[0]
        LocalDict={'collectionMembers':[],'references':{'name':[]}}
        Collection=False
        for SmallMolecule_property in single_SmallMolecule:
            Collection=MetaParser_SecLoop(LocalDict, SmallMolecule_property, Collection)
            if '}entityReference' in SmallMolecule_property.tag:
                LocalDict['references']=zipDicts(SmallMoleculeRefs[SmallMolecule_property.attrib.values()[0][1:]], LocalDict['references'])
        if Collection:
            SmallMolecule_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            SmallMolecules[key]=LocalDict

def parse_Proteins(root):
    Proteins_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Protein')
    for single_Protein in Proteins_Xml:
        key=single_Protein.attrib.values()[0]
        LocalDict={'collectionMembers':[],'modification':[],'references':{'name':[]}}
        Collection=False
        for Protein_property in single_Protein:
            Collection=MetaParser_SecLoop(LocalDict, Protein_property, Collection)
            if '}entityReference' in Protein_property.tag:
                LocalDict['references']=zipDicts(ProteinRefs[Protein_property.attrib.values()[0][1:]], LocalDict['references'])
            if '}feature' in Protein_property.tag and 'ModificationFeature' in Protein_property.attrib.values()[0]:
                LocalDict['modification'].append(ModificationFeatures[Protein_property.attrib.values()[0][1:]])
        if LocalDict['modification']==[]:
            del LocalDict['modification']
        if Collection:
            Protein_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            Proteins[key]=LocalDict
    
def parse_PhysicalEntities(root):
    PhysicalEntities_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity')
    for single_PhysicalEntity in PhysicalEntities_Xml:
        key=single_PhysicalEntity.attrib.values()[0]
        LocalDict={'collectionMembers':[],'references':{'name':[]}}
        Collection=False
        for PhysicalEntity_property in single_PhysicalEntity:
            Collection=MetaParser_SecLoop(LocalDict, PhysicalEntity_property, Collection)
        if Collection:
            PhysicalEntity_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            PhysicalEntities[key]=LocalDict
    
def parse_Complexes(root):
    Complex_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Complex')
    for single_Complex in Complex_Xml:
        key=single_Complex.attrib.values()[0]
        LocalDict={'collectionMembers':[],'parts':[],'references':{'name':[]}}
        Collection=False
        for complex_property in single_Complex:
            Collection=MetaParser_SecLoop(LocalDict, complex_property, Collection)
            if '}component' in complex_property.tag and not 'Stoichiometry' in complex_property.tag:
                LocalDict['parts'].append(complex_property.attrib.values()[0][1:])
        if Collection:
            del LocalDict['parts']
            Complex_Collections[key]=LocalDict
        else:
            del LocalDict['collectionMembers']
            Complexes[key]=LocalDict
    
def parse_TemplateReactions(root):
    TemplateReaction_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReaction')
    for single_TemplateReaction in TemplateReaction_Xml:
        key=single_TemplateReaction.attrib.values()[0]
        LocalDict={'references':{'name':[]}}
        for TemplateReaction_property in single_TemplateReaction:
            if '}product' in TemplateReaction_property.tag:
                LocalDict['product']=TemplateReaction_property.attrib.values()[0][1:]
            if '}displayName' in TemplateReaction_property.tag:
                LocalDict['displayName']=TemplateReaction_property.text
            if '}name' in TemplateReaction_property.tag:
                LocalDict['references']['name']=TemplateReaction_property.text
        TemplateReactions[key]=LocalDict
    
def parse_Degradations(root):
    DegradationXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Degradation')
    for single_Degradation in DegradationXml:
        key=single_Degradation.attrib.values()[0]
        LocalDict={'references':{'name':[]}}
        for Degradation_property in single_Degradation:
            if '}left' in Degradation_property.tag:
                LocalDict['product']=Degradation_property.attrib.values()[0][1:]
            if '}displayName' in Degradation_property.tag:
                LocalDict['displayName']=Degradation_property.text
            if '}eCNumber' in Degradation_property.tag:
                LocalDict['references']['eCNumber']=Degradation_property.text
        Degradations[key]=LocalDict
    
def parse_BiochemicalReactions(root):
    BiochemicalReactionXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}BiochemicalReaction')
    for single_BiochemicalReaction in BiochemicalReactionXml:
        key=single_BiochemicalReaction.attrib.values()[0]
        LocalDict={'right':[],'left':[],'references':{'name':[]}}
        for BiochemicalReaction_property in single_BiochemicalReaction:
            if '}left' in BiochemicalReaction_property.tag:
                LocalDict['left'].append(BiochemicalReaction_property.attrib.values()[0][1:])
            if '}right' in BiochemicalReaction_property.tag:
                LocalDict['right'].append(BiochemicalReaction_property.attrib.values()[0][1:])
            if '}displayName' in BiochemicalReaction_property.tag:
                LocalDict['displayName']=BiochemicalReaction_property.text
            if '}eCNumber' in BiochemicalReaction_property.tag:
                LocalDict['references']['eCNumber']=BiochemicalReaction_property.text
        BiochemicalReactions[key]=LocalDict
    
def parse_Catalysises(root):
    CatalysisesXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Catalysis')
    for single_Catalysis in CatalysisesXml:
        key=single_Catalysis.attrib.values()[0]
        Catalysises[key]={}
        for catalysis_property in single_Catalysis:
            if '}controlled' in catalysis_property.tag :
                Catalysises[key]['controlled']=catalysis_property.attrib.values()[0][1:]
            if '}controller' in catalysis_property.tag:
                Catalysises[key]['controller']=catalysis_property.attrib.values()[0][1:]
            if '}controlType' in catalysis_property.tag:
                Catalysises[key]['ControlType']=catalysis_property.text
    TemplateReactRegulationXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
    for single_TRR in TemplateReactRegulationXml:
        key=single_TRR.attrib.values()[0]
        Catalysises[key]={}
        for TRR_property in single_TRR:
            if '}controlled' in TRR_property.tag:
                Catalysises[key]['controlled']=TRR_property.attrib.values()[0][1:]
            if '}controller' in TRR_property.tag:
                Catalysises[key]['controller']=TRR_property.attrib.values()[0][1:]
            if '}controlType' in TRR_property.tag:
                Catalysises[key]['ControlType']=TRR_property.text
    ControlXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Control')
    for single_Control in ControlXml:
        key=single_Control.attrib.values()[0]
        Catalysises[key]={}
        for control_property in single_Control:
            if '}controlled' in control_property.tag and 'Pathway' not in control_property.attrib.values()[0][1:]:
                Catalysises[key]['controlled']=control_property.attrib.values()[0][1:]
            if '}controller' in control_property.tag:
                Catalysises[key]['controller']=control_property.attrib.values()[0][1:]
            if '}controlType' in control_property.tag and 'BiochemicalReaction' in control_property.attrib.values()[0]:
                Catalysises[key]['ControlType']=control_property.text
    
def parse_Modulations(root):
    ModulationsXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Modulation')
    for single_Modulation in ModulationsXml:
        key=single_Modulation.attrib.values()[0]
        Modulations[key]={}
        for modulation_property in single_Modulation:
            if '}displayName' in modulation_property.tag:
                Modulations[key]['displayName']=modulation_property.text
            if '}controller' in modulation_property.tag:
                Modulations[key]['controller']=modulation_property.attrib.values()[0][1:]
            if '}controlled' in modulation_property.tag:
                Modulations[key]['controlled']=Catalysises[modulation_property.attrib.values()[0][1:]]['controller']
            if '}controlType' in modulation_property.tag:
                Modulations[key]['controlType']=modulation_property.text

def parse_Pathways(root):
    PathwaysXml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}Pathway')
    for single_Pathway in PathwaysXml:
        key=single_Pathway.attrib.values()[0]
        LocalDict={'components':[],'references':{'name':[]}, 'PathwayStep':[]}
        for pathway_property in single_Pathway:
            if '}displayName' in pathway_property.tag:
                LocalDict['displayName']=pathway_property.text
            if '}name' in pathway_property.tag:
                LocalDict['references']['name'].append(pathway_property.text)
            if '}pathwayComponent' in pathway_property.tag and 'Pathway' in pathway_property.attrib.values()[0]:
                LocalDict['components'].append(pathway_property.attrib.values()[0][1:])
            if '}pathwayOrder' in pathway_property.tag:
                LocalDict['PathwayStep'].append(pathway_property.attrib.values()[0][1:])    
        Pathways[key]=LocalDict

def parse_Pathway_Steps(root):
    Pathway_Steps_Xml=root.findall('{http://www.biopax.org/release/biopax-level3.owl#}PathwayStep')
    exclude=['Modulation','Control','TemplateReactionRegulation', 'Catalysis']
    for single_Pathway_step in Pathway_Steps_Xml:
        key=single_Pathway_step.attrib.values()[0]
        LocalDict={'components':[], 'nextStep':[]}
        for pathway_property in single_Pathway_step:
            if '}stepProcess ' in pathway_property.tag and not any(x in pathway_property.attrib.values()[0] for x in exclude):
                LocalDict['components'].append(pathway_property.attrib.values()[0][1:])
            if '}nextStep ' in pathway_property.tag:
                LocalDict['nextStep'].append(pathway_property.attrib.values()[0][1:])    
        PathwaySteps[key]=LocalDict

def parse_all():


    ####################################################################################
    #
    # Initial parse of all the elements in the Reactome BioPax package
    #
    ####################################################################################

    tree = ET.parse(conf.ReactomeBioPax)
    root = tree.getroot()


    parse_BioSource(root)
    parse_CellularLocations(root)
    parse_SeqModVoc(root)
    parse_SeqSite(root)
    
    parse_DnaRefs(root)
    parse_RnaRefs(root)
    parse_SmallMoleculeRefs(root)
    parse_ProteinRefs(root)
    parse_ModificationFeatures(root)
    
    parse_Dnas(root)
    parse_Rnas(root)
    parse_SmallMolecules(root)
    parse_Proteins(root)
    parse_PhysicalEntities(root)
    parse_Complexes(root)
    
    parse_TemplateReactions(root)
    parse_Degradations(root)
    parse_BiochemicalReactions(root)
    parse_Catalysises(root)
    
    parse_Modulations(root)
    parse_Pathways(root)
    parse_Pathway_Steps(root)

parse_all()

if __name__ == "__main__":
    acnums = []
    for i, (protein, chardict) in enumerate(Proteins.iteritems()):
        if 'UniProt' in chardict['references'].keys():
            print i, chardict['references']['UniProt']
            acnums.append(chardict['references']['UniProt'])
    print list(set(acnums))
    print len(list(set(acnums)))