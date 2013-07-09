'''
Created on Jun 15, 2013

@author: andrei
'''
import logging
import Reactome_org_parser as DG
from neo4j_Declarations.Graph_Declarator import DatabaseGraph

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



LocalDict={} # accelerated access pointer to the objects

def InsertCellLocations():
    for Loc in DG.CellularLocations.keys():
        LocalDict[Loc]=DatabaseGraph.Location.create(ID=Loc, displayName=DG.CellularLocations[Loc])

def MinimalAnnotInsert(primary,reflist):
    for Type in reflist.keys():
        if Type!='name' and reflist[Type]!='' and reflist[Type]!=[]:
            secondary=DatabaseGraph.AnnotNode.create(ptype=Type,payload=reflist[Type])
            print primary, secondary, primary.ID
            DatabaseGraph.is_annotated.create(primary, secondary, costum_from=primary.ID,costum_to='Annotation')

def MetaInsert(function,dico):
    length=len(dico)
    counter=0
    for key in dico.keys():
        counter+=1
        print '\n',counter,'/',length,'\n'
        primary=function.create(ID=key, displayName=dico[key]['displayName'], localization=dico[key]['cellularLocation'])
        LocalDict[key]=primary
        MinimalAnnotInsert(LocalDict[key], dico[key]['references'])
        if 'cellularLocation' in dico[key].keys():
            secondary=LocalDict[dico[key]['cellularLocation']]
            DatabaseGraph.is_localized.create(primary, secondary,costum_from=primary.ID,costum_to=secondary.ID)
            # TODO add ModificationFeature insertion
        if 'modification' in dico[key].keys():
            for modification in dico[key]['modification']:
                if 'location' in modification.keys() and 'modification' in modification.keys():
                    LocMod=DatabaseGraph.ModificationFeature.create(ID = modification['ID'], type="post-translational_Mod", location=modification['location'], displayName=modification['modification'])
                    DatabaseGraph.is_able_to_modify.create(primary,LocMod,costum_from=primary.ID,costum_to=LocMod.ID)


def CollectionRefsInsert(primaryCollection):
    for key in primaryCollection.keys():
        for ref in primaryCollection[key]['collectionMembers']:
            DatabaseGraph.is_part_of_collection.create(LocalDict[key],LocalDict[ref],costum_from=LocalDict[key].ID,costum_to=LocalDict[ref].ID)

def ComplexPartsInsert():
    for key in DG.Complexes.keys():
        for part in DG.Complexes[key]['parts']:
            if 'Stoichiometry' not in part:
                DatabaseGraph.is_part_of_complex.create(LocalDict[key],LocalDict[part],costum_from=LocalDict[key].ID,costum_to=LocalDict[part].ID)

def ReactionInsert(function,dico):
    for key in dico.keys():
        LocalDict[key]=function.create(ID=key,displayName=dico[key]['displayName'])
        MinimalAnnotInsert(LocalDict[key], dico[key]['references'])
        for subkey in dico[key].keys():
            if subkey in ['left','right']:
                for elt in dico[key][subkey]:
                    DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[elt],side=subkey,costum_from=LocalDict[key].ID,costum_to=LocalDict[elt].ID) 
            if subkey=='product':
                DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[dico[key][subkey]],costum_from=LocalDict[key].ID,costum_to=LocalDict[dico[key][subkey]].ID)

def CatalysisInsert():
    for key in DG.Catalysises.keys():
        if 'controller' in DG.Catalysises[key].keys() and 'controlled' in DG.Catalysises[key].keys():
            if DG.Catalysises[key]['controlled'] in LocalDict.keys() and DG.Catalysises[key]['controller'] in LocalDict.keys():
                if 'ControlType' in DG.Catalysises[key].keys():
                    primary=LocalDict[DG.Catalysises[key]['controller']]
                    secondary=LocalDict[DG.Catalysises[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary, secondary, ID=key, controlType=DG.Catalysises[key]['ControlType'],costum_from=primary.ID,costum_to=secondary.ID)
                else:
                    primary=LocalDict[DG.Catalysises[key]['controller']]
                    secondary=LocalDict[DG.Catalysises[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary, secondary, ID=key, controlType='UNKNOWN',costum_from=primary.ID,costum_to=secondary.ID)
            else: 
                logging.debug("\t%s : %s, %s, %s", key, DG.Catalysises[key],DG.Catalysises[key]['controlled'] in LocalDict.keys(),DG.Catalysises[key]['controller'] in LocalDict.keys())
        else: 
            logging.debug("%s : %s, %s, %s,", key, DG.Catalysises[key],'controller' in DG.Catalysises[key].keys(),'controlled' in DG.Catalysises[key].keys())

def ModulationInsert():
    for key in DG.Modulations.keys():
        primary=LocalDict[DG.Modulations[key]['controller']]
        secondary=LocalDict[DG.Modulations[key]['controlled']]
        LocalDict[key]=DatabaseGraph.is_regulant.create(primary, secondary, ID=key,controlType=DG.Modulations[key]['controlType'],costum_from=primary.ID,costum_to=secondary.ID)

def Pathways_Insert():
    for key in DG.PathwaySteps.keys():
        primary=DatabaseGraph.PathwayStep.create(ID=key)
        LocalDict[key]=primary
    for key in DG.Pathways.keys():
        primary=DatabaseGraph.Pathway.create(ID=key, displayName=DG.Pathways['displayName'])
        LocalDict[key]=primary
    for key in DG.PathwaySteps.keys():
        for component in DG.PathwaySteps[key]['components']:
            primary=LocalDict[key]
            secondary=LocalDict[component]
            DatabaseGraph.is_part_of_pathway.create(primary,secondary,costum_from=key,costum_to=component)
        for nextStep in DG.PathwaySteps[key]['nextStep']:
            primary=LocalDict[key]
            secondary=LocalDict[nextStep]
            DatabaseGraph.is_next_in_pathway.create(primary,secondary,costum_from=key,costum_to=nextStep)
    for key in DG.Pathways.keys():
        for pathwayStep in DG.Pathways[key]['PathwayStep']:
            primary=LocalDict[key]
            secondary=LocalDict[pathwayStep]
            DatabaseGraph.is_part_of_pathway.create(primary,secondary,costum_from=key,costum_to=pathwayStep)
        for Sub_Pathway in DG.Pathways[key]['components']:
            primary=LocalDict[key]
            secondary=LocalDict[Sub_Pathway]
            DatabaseGraph.is_part_of_pathway.create(primary,secondary,costum_from=key,costum_to=Sub_Pathway)

def getOneMetaSet(function):
    for MetaKey in function.get_all():
        if MetaKey!=None:
            LocalDict[MetaKey.ID]=MetaKey
        
def getAllMetaSets():
    functionList=[
                  DatabaseGraph.DNA,
                  DatabaseGraph.DNA_Collection,
                  DatabaseGraph.RNA,
                  DatabaseGraph.RNA_Collection,
                  DatabaseGraph.SmallMolecule,
                  DatabaseGraph.SmallMolecule_Collection,
                  DatabaseGraph.Protein,
                  DatabaseGraph.Protein_Collection,
                  DatabaseGraph.Complex,
                  DatabaseGraph.Complex_Collection,
                  DatabaseGraph.PhysicalEntity,
                  DatabaseGraph.PhysicalEntity_Collection
                  ]
    
    for function in functionList:
        getOneMetaSet(function)

#
# InsertCellLocations()
# 
# MetaInsert(DatabaseGraph.DNA, DG.Dnas)
# MetaInsert(DatabaseGraph.DNA_Collection, DG.Dna_Collections)
# MetaInsert(DatabaseGraph.RNA, DG.Rnas)
# MetaInsert(DatabaseGraph.RNA_Collection, DG.Rna_Collections)
# MetaInsert(DatabaseGraph.SmallMolecule, DG.SmallMolecules)
# MetaInsert(DatabaseGraph.SmallMolecule_Collection, DG.SmallMolecule_Collections)
# MetaInsert(DatabaseGraph.Protein, DG.Proteins)
# MetaInsert(DatabaseGraph.Protein_Collection, DG.Protein_Collections)
# MetaInsert(DatabaseGraph.Complex, DG.Complexes)
# MetaInsert(DatabaseGraph.Complex_Collection, DG.Complex_Collections)
# MetaInsert(DatabaseGraph.PhysicalEntity, DG.PhysicalEntities)
# MetaInsert(DatabaseGraph.PhysicalEntity_Collection, DG.PhysicalEntity_Collections)
# CollectionRefsInsert(DG.Dna_Collections)
# CollectionRefsInsert(DG.Rna_Collections)
# CollectionRefsInsert(DG.SmallMolecule_Collections)
# CollectionRefsInsert(DG.Protein_Collections)
# CollectionRefsInsert(DG.Complex_Collections)
# CollectionRefsInsert(DG.PhysicalEntity_Collections)
#
# ComplexPartsInsert()
# 
# ## Meta insert finished
# ReactionInsert(DatabaseGraph.TemplateReaction, DG.TemplateReactions)
# ReactionInsert(DatabaseGraph.Degradation, DG.Degradations)
# ReactionInsert(DatabaseGraph.BiochemicalReaction, DG.BiochemicalReactions)
# 
# ## Reaction insert finished
# CatalysisInsert()
# ModulationInsert()
# Pathways_Insert()
#

