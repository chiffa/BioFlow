'''
Created on Jun 15, 2013

@author: andrei
'''
import logging
import Reactome_org_parser as DG
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph

####################################################################################
#
# Logger behavior definition (most of the time it fails to function due to the)
# logs collision with neo4j
#
####################################################################################

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

def MinimalAnnotInsert(annotated_node, payload_list):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct, local DB ID and thus)
    needs to be inserted at the same time as the annotated object
    """
    for Type in payload_list.keys():
        if Type!='name' and payload_list[Type]!='' and payload_list[Type]!=[]:
            if type(Type)!=list:
                secondary=DatabaseGraph.AnnotNode.create(ptype=Type,payload=payload_list[Type])
                print annotated_node, secondary, annotated_node.ID
                DatabaseGraph.is_annotated.create(annotated_node, secondary, costum_from=annotated_node.ID,costum_to='Annotation')
            else:
                for subelt in payload_list[Type]:
                    secondary=DatabaseGraph.AnnotNode.create(ptype=Type,payload=subelt)
                    print annotated_node, secondary, annotated_node.ID
                    DatabaseGraph.is_annotated.create(annotated_node, secondary, costum_from=annotated_node.ID,costum_to='Annotation')

def MetaInsert(bulbs_graph_class, property_source_dict):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a bulbs class and pumping the
    source information from the property source
    """
    length=len(property_source_dict)
    counter=0
    for key in property_source_dict.keys():
        counter+=1
        print '\n',counter,'/',length,'\n'
        primary=bulbs_graph_class.create(ID=key, displayName=property_source_dict[key]['displayName'], localization=property_source_dict[key]['cellularLocation'])
        LocalDict[key]=primary
        MinimalAnnotInsert(LocalDict[key], property_source_dict[key]['references'])
        if 'cellularLocation' in property_source_dict[key].keys():
            secondary=LocalDict[property_source_dict[key]['cellularLocation']]
            DatabaseGraph.is_localized.create(primary, secondary,costum_from=primary.ID,costum_to=secondary.ID)
            # TODO add ModificationFeature insertion
        if 'modification' in property_source_dict[key].keys():
            for modification in property_source_dict[key]['modification']:
                if 'location' in modification.keys() and 'modification' in modification.keys():
                    LocMod=DatabaseGraph.ModificationFeature.create(ID = modification['ID'], type="post-translational_Mod", location=modification['location'], displayName=modification['modification'])
                    DatabaseGraph.is_able_to_modify.create(primary,LocMod,costum_from=primary.ID,costum_to=LocMod.ID)

def CollectionRefsInsert(primaryCollection):
    """
    Links a collection object reference to the members of the collection.
    """
    for key in primaryCollection.keys():
        for ref in primaryCollection[key]['collectionMembers']:
            DatabaseGraph.is_part_of_collection.create(LocalDict[key],LocalDict[ref],costum_from=LocalDict[key].ID,costum_to=LocalDict[ref].ID)

def ComplexPartsInsert():
    """
    Links part of a complex to the complex
    """
    for key in DG.Complexes.keys():
        for part in DG.Complexes[key]['parts']:
            if 'Stoichiometry' not in part:
                DatabaseGraph.is_part_of_complex.create(LocalDict[key],LocalDict[part],costum_from=LocalDict[key].ID,costum_to=LocalDict[part].ID)

def ReactionInsert(bulbs_graph_class, property_source_dict):
    """
    Inserst a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs class and pumping the
    source information from the property source
    """
    for key in property_source_dict.keys():
        LocalDict[key]=bulbs_graph_class.create(ID=key,displayName=property_source_dict[key]['displayName'])
        MinimalAnnotInsert(LocalDict[key], property_source_dict[key]['references'])
        for subkey in property_source_dict[key].keys():
            if subkey in ['left','right']:
                for elt in property_source_dict[key][subkey]:
                    DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[elt],side=subkey,costum_from=LocalDict[key].ID,costum_to=LocalDict[elt].ID) 
            if subkey=='product':
                DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[property_source_dict[key][subkey]],costum_from=LocalDict[key].ID,costum_to=LocalDict[property_source_dict[key][subkey]].ID)

def CatalysisInsert():
    """
    Inserts all the catalysis links from one meta-element to an another
    """
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
    """
    Inserts all the Modulation links from one meta-element to an another
    """
    for key in DG.Modulations.keys():
        primary=LocalDict[DG.Modulations[key]['controller']]
        secondary=LocalDict[DG.Modulations[key]['controlled']]
        LocalDict[key]=DatabaseGraph.is_regulant.create(primary, secondary, ID=key,controlType=DG.Modulations[key]['controlType'],costum_from=primary.ID,costum_to=secondary.ID)

def Pathways_Insert():
    """
    Inserts all the Pathways, linking and chaining subpathways
    Attention, it have to be imported at the same time as the reactions.
    """
    for key in DG.PathwaySteps.keys():
        primary=DatabaseGraph.PathwayStep.create(ID=key)
        LocalDict[key]=primary
    for key in DG.Pathways.keys():
        primary=DatabaseGraph.Pathway.create(ID=key, displayName=DG.Pathways[key]['displayName'])
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
    """
    In case a MetaObject was already inserted, reloads it to the local dictionary for futher annoation
    """
    for MetaKey in function.get_all():
        if MetaKey!=None:
            LocalDict[MetaKey.ID]=MetaKey
        
def getAllMetaSets():
    """
    In case the MetaObjects were already inserted, reloads them all to the local dictionary for futher annoation
    """
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

# Attention, thuis function has to be redifined regularly, since it is mediated by the
def clean():
    """
    Serves to clear a specific set and nodes in case an insertion went awry
    """
    for PathStep in DatabaseGraph.PathwayStep.get_all():
        ID=str(PathStep).split('/')[-1][:-1]
        DatabaseGraph.PathwayStep.delete(ID)


def insert(Skip):
    """
    Performs the massive import of the Reactome database into the local neo4j database.
    :parameter Skip:
    * N => will skip nothing and implement the import once and for all.
    * M => skips meta import, recovers the metas and resumes from the Reactions import.
    """
    # TODO: implement the skipping behavior

    if Skip=='N':
        InsertCellLocations()

        MetaInsert(DatabaseGraph.DNA, DG.Dnas)
        MetaInsert(DatabaseGraph.DNA_Collection, DG.Dna_Collections)
        MetaInsert(DatabaseGraph.RNA, DG.Rnas)
        MetaInsert(DatabaseGraph.RNA_Collection, DG.Rna_Collections)
        MetaInsert(DatabaseGraph.SmallMolecule, DG.SmallMolecules)
        MetaInsert(DatabaseGraph.SmallMolecule_Collection, DG.SmallMolecule_Collections)
        MetaInsert(DatabaseGraph.Protein, DG.Proteins)
        MetaInsert(DatabaseGraph.Protein_Collection, DG.Protein_Collections)
        MetaInsert(DatabaseGraph.Complex, DG.Complexes)
        MetaInsert(DatabaseGraph.Complex_Collection, DG.Complex_Collections)
        MetaInsert(DatabaseGraph.PhysicalEntity, DG.PhysicalEntities)
        MetaInsert(DatabaseGraph.PhysicalEntity_Collection, DG.PhysicalEntity_Collections)
        CollectionRefsInsert(DG.Dna_Collections)
        CollectionRefsInsert(DG.Rna_Collections)
        CollectionRefsInsert(DG.SmallMolecule_Collections)
        CollectionRefsInsert(DG.Protein_Collections)
        CollectionRefsInsert(DG.Complex_Collections)
        CollectionRefsInsert(DG.PhysicalEntity_Collections)

        ComplexPartsInsert()



    if Skip == 'M':
        getAllMetaSets()

    ## Meta insert/retrieval finished

    ReactionInsert(DatabaseGraph.TemplateReaction, DG.TemplateReactions)
    ReactionInsert(DatabaseGraph.Degradation, DG.Degradations)
    ReactionInsert(DatabaseGraph.BiochemicalReaction, DG.BiochemicalReactions)

    ## Reaction insert finished
    CatalysisInsert()
    ModulationInsert()
    Pathways_Insert()

if __name__ == "__main__":
    insert('N')