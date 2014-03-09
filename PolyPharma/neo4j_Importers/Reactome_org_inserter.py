'''
Created on Jun 15, 2013

@author: andrei
'''
import logging
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
import pickle
from PolyPharma.configs import Dumps, Leg_ID_Filter

####################################################################################
#
# Logger behavior definition (most of the time it fails to function due to the)
# logs collision with neo4j
#
####################################################################################

# TODO: export logs location to the configs file
# TODO: create a dict that is pickled into ther reserve to remember the forbidden IDs

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='../logs/dynamics_full.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

LocalDict = {} # accelerated access pointer to the objects
ForbiddenIDs = []


def InsertCellLocations(Cell_locations_dict):
    for Loc in Cell_locations_dict.keys():
        LocalDict[Loc] = DatabaseGraph.Location.create(ID = Loc,
                                                     displayName = Cell_locations_dict[Loc])

def MinimalAnnotInsert(annotated_node, payload_list):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct, local DB ID and thus)
    needs to be inserted at the same time as the annotated object
    """
    for Type in payload_list.keys():
        if Type!='name' and payload_list[Type]!='' and payload_list[Type]!=[]:
            if type(Type)!=list:
                secondary=DatabaseGraph.AnnotNode.create(ptype=Type,
                                                         payload=payload_list[Type])
                print annotated_node, secondary, annotated_node.ID
                DatabaseGraph.is_annotated.create(annotated_node,
                                                  secondary,
                                                  costum_from=annotated_node.ID,
                                                  costum_to='Annotation')
            else:
                for subelt in payload_list[Type]:
                    secondary=DatabaseGraph.AnnotNode.create(ptype=Type,payload=subelt)
                    print annotated_node, secondary, annotated_node.ID
                    DatabaseGraph.is_annotated.create(annotated_node,
                                                      secondary,
                                                      costum_from=annotated_node.ID,
                                                      costum_to='Annotation')

def MetaInsert(bulbs_graph_class, property_source_dict):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a bulbs class and pumping the
    source information from the property source
    """
    length = len(property_source_dict)
    counter = 0
    for key in property_source_dict.keys():
        counter += 1
        print '\n', counter, '/', length, '\n'
        primary = bulbs_graph_class.create(ID = key,
                                         displayName = property_source_dict[key]['displayName'],
                                         localization = property_source_dict[key]['cellularLocation'],
                                         main_connex = False)
        if key in Leg_ID_Filter:
            ForbiddenIDs.append(str(primary).split('/')[-1][:-1])
        LocalDict[key] = primary
        MinimalAnnotInsert(LocalDict[key], property_source_dict[key]['references'])
        if 'cellularLocation' in property_source_dict[key].keys():
            secondary = LocalDict[property_source_dict[key]['cellularLocation']]
            DatabaseGraph.is_localized.create(primary,
                                              secondary,
                                              costum_from = primary.ID,
                                              costum_to = secondary.ID)
        if 'modification' in property_source_dict[key].keys():
            for modification in property_source_dict[key]['modification']:
                if 'location' in modification.keys() and 'modification' in modification.keys():
                    LocMod=DatabaseGraph.ModificationFeature.create(ID = modification['ID'],
                                                                    type = "post-translational_Mod",
                                                                    location = modification['location'],
                                                                    displayName = modification['modification'])
                    DatabaseGraph.is_able_to_modify.create(primary,
                                                           LocMod,
                                                           costum_from = primary.ID,
                                                           costum_to = LocMod.ID)

def CollectionRefsInsert(primaryCollection):
    """
    Links a collection object reference to the members of the collection.
    """
    for key in primaryCollection.keys():
        for ref in primaryCollection[key]['collectionMembers']:
            DatabaseGraph.is_part_of_collection.create(LocalDict[key],
                                                       LocalDict[ref],
                                                       costum_from=LocalDict[key].ID,
                                                       costum_to=LocalDict[ref].ID)

def ComplexPartsInsert(complexes_dict):
    """
    Links part of a complex to the complex
    """
    for key in complexes_dict.keys():
        for part in complexes_dict[key]['parts']:
            if 'Stoichiometry' not in part:
                DatabaseGraph.is_part_of_complex.create(LocalDict[key],
                                                        LocalDict[part],
                                                        costum_from=LocalDict[key].ID,
                                                        costum_to=LocalDict[part].ID)

def ReactionInsert(bulbs_graph_class, property_source_dict):
    """
    Inserst a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs class and pumping the
    source information from the property source
    """
    for key in property_source_dict.keys():
        LocalDict[key]=bulbs_graph_class.create(ID=key,
                                                displayName=property_source_dict[key]['displayName'])
        MinimalAnnotInsert(LocalDict[key], property_source_dict[key]['references'])
        for subkey in property_source_dict[key].keys():
            if subkey in ['left','right']:
                for elt in property_source_dict[key][subkey]:
                    DatabaseGraph.is_reaction_participant.create(LocalDict[key],
                                                                 LocalDict[elt],
                                                                 side = subkey,
                                                                 costum_from = LocalDict[key].ID,
                                                                 costum_to = LocalDict[elt].ID)
            if subkey=='product':
                DatabaseGraph.is_reaction_participant.create(LocalDict[key],
                                                             LocalDict[property_source_dict[key][subkey]],
                                                             costum_from = LocalDict[key].ID,
                                                             costum_to = LocalDict[property_source_dict[key][subkey]].ID)

def CatalysisInsert(catalysises_dict):
    """
    Inserts all the catalysis links from one meta-element to an another
    """
    for key in catalysises_dict.keys():
        if 'controller' in catalysises_dict[key].keys() and 'controlled' in catalysises_dict[key].keys():
            if catalysises_dict[key]['controlled'] in LocalDict.keys() and catalysises_dict[key]['controller'] in LocalDict.keys():
                if 'ControlType' in catalysises_dict[key].keys():
                    primary=LocalDict[catalysises_dict[key]['controller']]
                    secondary=LocalDict[catalysises_dict[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary,
                                                                      secondary,
                                                                      ID=key,
                                                                      controlType=catalysises_dict[key]['ControlType'],
                                                                      costum_from=primary.ID,
                                                                      costum_to=secondary.ID)
                else:
                    primary=LocalDict[catalysises_dict[key]['controller']]
                    secondary=LocalDict[catalysises_dict[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary,
                                                                      secondary,
                                                                      ID=key,
                                                                      controlType='UNKNOWN',
                                                                      costum_from=primary.ID,
                                                                      costum_to=secondary.ID)
            else: 
                logging.debug("\t%s : %s, %s, %s", key, catalysises_dict[key],
                              catalysises_dict[key]['controlled'] in LocalDict.keys(),
                              catalysises_dict[key]['controller'] in LocalDict.keys())
        else: 
            logging.debug("%s : %s, %s, %s,", key, catalysises_dict[key],
                          'controller' in catalysises_dict[key].keys(),
                          'controlled' in catalysises_dict[key].keys())

def ModulationInsert(modulations_dict):
    """
    Inserts all the Modulation links from one meta-element to an another
    """
    for key in modulations_dict.keys():
        primary=LocalDict[modulations_dict[key]['controller']]
        secondary=LocalDict[modulations_dict[key]['controlled']]
        LocalDict[key]=DatabaseGraph.is_regulant.create(primary,
                                                        secondary,
                                                        ID = key,
                                                        controlType = modulations_dict[key]['controlType'],
                                                        costum_from = primary.ID,
                                                        costum_to = secondary.ID)

def Pathways_Insert(pathway_steps, pathways):
    """
    Inserts all the Pathways, linking and chaining subpathways
    Attention, it have to be imported at the same time as the reactions.
    """
    for key in pathway_steps.keys():
        primary=DatabaseGraph.PathwayStep.create(ID=key)
        LocalDict[key]=primary
    for key in pathways.keys():
        primary=DatabaseGraph.Pathway.create(ID=key,
                                             displayName=pathways[key]['displayName'])
        LocalDict[key]=primary
    for key in pathway_steps.keys():
        for component in pathway_steps[key]['components']:
            primary=LocalDict[key]
            secondary=LocalDict[component]
            DatabaseGraph.is_part_of_pathway.create(primary,
                                                    secondary,
                                                    costum_from=key,
                                                    costum_to=component)
        for nextStep in pathway_steps[key]['nextStep']:
            primary=LocalDict[key]
            secondary=LocalDict[nextStep]
            DatabaseGraph.is_next_in_pathway.create(primary,
                                                    secondary,
                                                    costum_from=key,
                                                    costum_to=nextStep)
    for key in pathways.keys():
        for pathwayStep in pathways[key]['PathwayStep']:
            primary=LocalDict[key]
            secondary=LocalDict[pathwayStep]
            DatabaseGraph.is_part_of_pathway.create(primary,
                                                    secondary,
                                                    costum_from=key,
                                                    costum_to=pathwayStep)
        for Sub_Pathway in pathways[key]['components']:
            primary=LocalDict[key]
            secondary=LocalDict[Sub_Pathway]
            DatabaseGraph.is_part_of_pathway.create(primary,
                                                    secondary,
                                                    costum_from=key,
                                                    costum_to=Sub_Pathway)


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
                  DatabaseGraph.PhysicalEntity_Collection,
                  ]
    
    for function in functionList:
        getOneMetaSet(function)


def clear_all(instruction_dict):
    """
    empties the whole BioPax-bound node set.
    """
    for name, bulbs_class in instruction_dict.iteritems():
        counter = 0
        if bulbs_class.get_all():
            for bulbs_class_instance in bulbs_class.get_all():
                counter += 1
                ID = str(bulbs_class_instance).split('/')[-1][:-1]
                DatabaseGraph.PathwayStep.delete(ID)
                if counter % 100 == 0:
                    print name, ':', counter

def run_diagnostics(instruction_dict):
    """
    Checks the number of nodes of each type.
    """
    supercounter = 0
    for name, bulbs_class in instruction_dict.iteritems():
        counter = 0
        if bulbs_class.get_all():
            for bulbs_class_instance in bulbs_class.get_all():
                counter += 1
        print name, ':', counter
        supercounter += counter
    print 'Total: ', supercounter


def insert_all(Skip='N'):
    """
    Performs the massive import of the Reactome database into the local neo4j database.
    :parameter Skip:
    * N => will skip nothing and implement the import once and for all.
    * M => skips meta import, recovers the metas and resumes from the Reactions import.
    """
    import Reactome_org_parser as DG

    if Skip=='N':
        InsertCellLocations(DG.CellularLocations)

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

        ComplexPartsInsert(DG.Complexes)

        # NOW dump the ForbiddenIDs
        pickle.dump(ForbiddenIDs, Dumps.Forbidden_IDs)

    if Skip == 'M':
        getAllMetaSets()

    ## Meta insert/retrieval finished

    ReactionInsert(DatabaseGraph.TemplateReaction, DG.TemplateReactions)
    ReactionInsert(DatabaseGraph.Degradation, DG.Degradations)
    ReactionInsert(DatabaseGraph.BiochemicalReaction, DG.BiochemicalReactions)

    ## Reaction insert finished
    CatalysisInsert(DG.Catalysises)
    ModulationInsert(DG.Modulations)
    Pathways_Insert(DG.PathwaySteps, DG.Pathways)



full_dict = {'DNA':DatabaseGraph.DNA,
              'DNA Collection':DatabaseGraph.DNA_Collection,
              'RNA':DatabaseGraph.RNA,
              'RNA Collection':DatabaseGraph.RNA_Collection,
              'Small Molecule':DatabaseGraph.SmallMolecule,
              'Small Molecule Collection':DatabaseGraph.SmallMolecule_Collection,
              'Protein':DatabaseGraph.Protein,
              'Protein Collection':DatabaseGraph.Protein_Collection,
              'Complex':DatabaseGraph.Complex,
              'Complex Collection':DatabaseGraph.Complex_Collection,
              'Physical Entity':DatabaseGraph.PhysicalEntity,
              'Physical Entity Collection':DatabaseGraph.PhysicalEntity_Collection,
              'TemplateReaction':DatabaseGraph.TemplateReaction,
              'Degradation':DatabaseGraph.Degradation,
              'BiochemicalReaction':DatabaseGraph.BiochemicalReaction,
              'Pathway Step':DatabaseGraph.PathwayStep,
              'Pathway':DatabaseGraph.Pathway,
              'Cell Locations':DatabaseGraph.Location,
              'Annotations':DatabaseGraph.AnnotNode,
              'Modification Feature':DatabaseGraph.ModificationFeature,
              'UNIPROT':DatabaseGraph.UNIPORT,
              'GO Term':DatabaseGraph.GOTerm,
            }

if __name__ == "__main__":
    # insert_all()
    # run_diagnostics(full_dict)
    # clear_all(full_dict)
    pass