"""
Created on Jun 15, 2013
:author: andrei
"""
import pickle
from BioFlow.utils.log_behavior import get_logger
from BioFlow.main_configs import Dumps, ReactomeBioPax
from BioFlow.internal_configs import Leg_ID_Filter
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph
from BioFlow.neo4j_db.db_io_routines import get_bulbs_id
from BioFlow.bio_db_parsers.reactomeParser import ReactomeParser

log = get_logger(__name__)

# TODO: refactor everything that uses memoization dictionary into a class
memoization_dict = {}  # accelerated access pointer to the objects
ForbiddenIDs = []


def insert_cell_locations(cell_locations_dict):
    """
    Creates nodes corresponding to cell locations

    :param cell_locations_dict:
    """
    for Loc in cell_locations_dict.keys():
        memoization_dict[Loc] = DatabaseGraph.Location.create(
            ID=Loc, displayName=cell_locations_dict[Loc])


def insert_minimal_annotations(annotated_node, annot_type_2_annot):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct,
    local DB ID and thus) needs to be inserted at the same time as the annotated object

    :param annotated_node:
    :param annot_type_2_annot:
    """
    for annotation_type, annotation_set in annot_type_2_annot.iteritems():
        if annotation_type != 'name' and annotation_set != '' and annotation_set != []:
            if not isinstance(annotation_type, list):
                annotation_set = [annotation_set]
            for annotation in annotation_set:
                annotation_node = DatabaseGraph.AnnotNode.create(
                    ptype=annotation_type, payload=annotation)
                log.info('created annotation %s for %s, _id:%s',
                         annotation_node, annotated_node, annotated_node.ID)
                DatabaseGraph.is_annotated.create(
                    annotated_node,
                    annotation_node,
                    costum_from=annotated_node.ID,
                    costum_to='Annotation')


def insert_meta_objects(bulbs_graph_class, meta_id_2_property_dict):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a
     bulbs class and pumping the source information from the property source

    :param bulbs_graph_class:
    :param meta_id_2_property_dict:
    """
    total_properties = len(meta_id_2_property_dict)

    for i, (meta_name, property_dict) in enumerate(meta_id_2_property_dict.iteritems()):
        log.info('meta_insert: property %s out of %s', i, total_properties)

        primary = bulbs_graph_class.create(
            ID=meta_name,
            displayName=property_dict['displayName'],
            localization=property_dict['cellularLocation'],
            main_connex=False)

        if meta_name in Leg_ID_Filter:
            ForbiddenIDs.append(get_bulbs_id(primary))

        memoization_dict[meta_name] = primary

        insert_minimal_annotations(
            memoization_dict[meta_name],
            meta_id_2_property_dict[meta_name]['references'])

        if 'cellularLocation' in meta_id_2_property_dict[meta_name].keys():
            secondary = memoization_dict[
                meta_id_2_property_dict[meta_name]['cellularLocation']]
            DatabaseGraph.is_localized.create(primary,
                                              secondary,
                                              costum_from=primary.ID,
                                              costum_to=secondary.ID)

        if 'modification' in meta_id_2_property_dict[meta_name].keys():
            for modification in meta_id_2_property_dict[meta_name]['modification']:
                if 'location' in modification.keys() and 'modification' in modification.keys():

                    located_modification = DatabaseGraph.ModificationFeature.create(
                        ID=modification['ID'],
                        type="post-translational_Mod",
                        location=modification['location'],
                        displayName=modification['modification'])

                    DatabaseGraph.is_able_to_modify.create(
                        primary, located_modification,
                        costum_from=primary.ID, costum_to=located_modification.ID)


def insert_collections(collections_2_members):
    """
    Links a collection object reference to the members of the collection.

    :param collections_2_members:
    """
    for collection, collection_property_dict in collections_2_members.iteritems():
        # TODO: add the display name
        for member in collection_property_dict['collectionMembers']:

            DatabaseGraph.is_part_of_collection.create(
                memoization_dict[collection],
                memoization_dict[member],
                costum_from=memoization_dict[collection].ID,
                costum_to=memoization_dict[member].ID)


def insert_complex_parts(complex_property_dict):
    """
    Links part of a complex to the complex

    :param complex_property_dict:
    """
    for key in complex_property_dict.keys():
        for part in complex_property_dict[key]['parts']:
            # TODO: remove redundant protection from Stoichiometry
            if 'Stoichiometry' not in part:
                DatabaseGraph.is_part_of_complex.create(
                    memoization_dict[key],
                    memoization_dict[part],
                    costum_from=memoization_dict[key].ID,
                    costum_to=memoization_dict[part].ID)


def insert_reactions(bulbs_graph_class, property_source_dict):
    """
    Inserts a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs
    class and pumping the source information from the property source

    :param bulbs_graph_class:
    :param property_source_dict:
    """
    for reaction, reaction_properties in property_source_dict.iteritems():
        memoization_dict[reaction] = bulbs_graph_class.create(
            ID=reaction, displayName=reaction_properties['displayName'])
        insert_minimal_annotations(
            memoization_dict[reaction],
            reaction_properties['references'])
        for property_name, property_value_list in reaction_properties.iteritems():
            if property_name in ['left', 'right']:
                for elt in property_value_list:
                    DatabaseGraph.is_reaction_participant.create(
                        memoization_dict[reaction],
                        memoization_dict[elt],
                        side=property_name,
                        costum_from=memoization_dict[reaction].ID,
                        costum_to=memoization_dict[elt].ID)


# TODO: catalysis and modulation are identical in the marked lines
def insert_catalysis(catalysises_dict):
    """
    Inserts all the catalysis links from one meta-element to an another

    :param catalysises_dict:
    """
    for catalysis, catalysis_properties in catalysises_dict.iteritems():

        if 'controller' in catalysis_properties.keys() \
                and 'controlled' in catalysis_properties.keys():

            if catalysis_properties['controlled'] in memoization_dict.keys() \
                    and catalysis_properties['controller'] in memoization_dict.keys():

                # TODO: this can be moved to the parsing
                if 'ControlType' not in catalysises_dict[catalysis].keys():
                    catalysis_properties['ControlType'] = 'UNKNOWN'

                controller = memoization_dict[catalysis_properties['controller']]  #
                controlled = memoization_dict[catalysis_properties['controlled']]  #
                memoization_dict[catalysis] = DatabaseGraph.is_catalysant.create(  #
                    controller,
                    controlled,
                    ID=catalysis,
                    controlType=catalysis_properties['ControlType'],
                    costum_from=controller.ID,
                    costum_to=controlled.ID)

            else:
                log.debug("Catalysis targets not memoized: %s : %s, %s, %s", catalysis,
                          catalysises_dict[catalysis],
                          catalysises_dict[catalysis]['controlled'] in memoization_dict.keys(),
                          catalysises_dict[catalysis]['controller'] in memoization_dict.keys())
        else:
            log.debug("Catalysis without control/controlled %s : %s, %s, %s,",
                      catalysis, catalysises_dict[catalysis],
                      'controller' in catalysises_dict[catalysis].keys(),
                      'controlled' in catalysises_dict[catalysis].keys())


def insert_modulation(modulations_dict):
    """
    Inserts all the Modulation links from one meta-element to an another

    :param modulations_dict:
    """
    for modulation, modulation_property_dict in modulations_dict.iteritems():
        controller = memoization_dict[modulation_property_dict['controller']]  #
        controlled = memoization_dict[modulation_property_dict['controlled']]  #
        memoization_dict[modulation] = DatabaseGraph.is_regulant.create(       #
            controller,
            controlled,
            ID=modulation,
            controlType=modulation_property_dict['controlType'],
            costum_from=controller.ID,
            costum_to=controlled.ID)


def insert_pathways(pathway_steps, pathways):
    """
    Inserts all the Pathways, linking and chaining subpathways
    Attention, it have to be imported at the same time as the reactions.

    :param pathway_steps:
    :param pathways:
    """
    for pathway_step in pathway_steps.keys():
        memoization_dict[pathway_step] = DatabaseGraph.PathwayStep.create(ID=pathway_step)

    for pathway_step in pathways.keys():
        memoization_dict[pathway_step] = DatabaseGraph.Pathway.create(
            ID=pathway_step, displayName=pathways[pathway_step]['displayName'])

    for pathway_step in pathway_steps.keys():

        for component in pathway_steps[pathway_step]['components']:
            DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                    memoization_dict[component],
                                                    costum_from=pathway_step,
                                                    costum_to=component)

        for next_step in pathway_steps[pathway_step]['nextStep']:
            DatabaseGraph.is_next_in_pathway.create(memoization_dict[pathway_step],
                                                    memoization_dict[next_step],
                                                    costum_from=pathway_step,
                                                    costum_to=next_step)
    for pathway_step in pathways.keys():

        for second_pathway_step in pathways[pathway_step]['PathwayStep']:
            DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                    memoization_dict[second_pathway_step],
                                                    costum_from=pathway_step,
                                                    costum_to=second_pathway_step)

        for sub_pathway in pathways[pathway_step]['components']:
            DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                    memoization_dict[sub_pathway],
                                                    costum_from=pathway_step,
                                                    costum_to=sub_pathway)


def get_one_meta_set(bulbs_bound_class):
    """
    In case a MetaObject was already inserted, reloads it to the local dictionary for further
    annotation insertion

    :param bulbs_bound_class:
    """
    for bulbs_object in bulbs_bound_class.get_all():
        if bulbs_object is not None:
            memoization_dict[bulbs_object.ID] = bulbs_object


def get_all_meta_sets():
    """
    In case the MetaObjects were already inserted, reloads them all to the local dictionary for
    further annotation insertion
    """
    # TODO: refactor using abbreviations and bulbs_names_dict
    list_of_bulbs_classes = [
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

    for function in list_of_bulbs_classes:
        get_one_meta_set(function)


def insert_reactome(skip_import='N'):
    """
    Performs the massive import of the Reactome database into the local neo4j database.

    :param skip_import:     * N => will skip nothing and implement the import once and for all.
                     * M => skips meta import, recovers the metas and resumes from the Reactions
                     import.
    """
    reactome_parser = ReactomeParser(ReactomeBioPax)
    reactome_parser.parse_all()

    # TODO: refactor using abbreviations and bulbs_names_dict
    if skip_import == 'N':

        insert_cell_locations(reactome_parser.CellularLocations)

        insert_meta_objects(DatabaseGraph.DNA, reactome_parser.Dnas)
        insert_meta_objects(DatabaseGraph.DNA_Collection, reactome_parser.Dna_Collections)
        insert_meta_objects(DatabaseGraph.RNA, reactome_parser.Rnas)
        insert_meta_objects(DatabaseGraph.RNA_Collection, reactome_parser.Rna_Collections)
        insert_meta_objects(DatabaseGraph.SmallMolecule, reactome_parser.SmallMolecules)
        insert_meta_objects(
            DatabaseGraph.SmallMolecule_Collection,
            reactome_parser.SmallMolecule_Collections)
        insert_meta_objects(DatabaseGraph.Protein, reactome_parser.Proteins)
        insert_meta_objects(DatabaseGraph.Protein_Collection, reactome_parser.Protein_Collections)
        insert_meta_objects(DatabaseGraph.Complex, reactome_parser.Complexes)
        insert_meta_objects(DatabaseGraph.Complex_Collection, reactome_parser.Complex_Collections)
        insert_meta_objects(DatabaseGraph.PhysicalEntity, reactome_parser.PhysicalEntities)
        insert_meta_objects(
            DatabaseGraph.PhysicalEntity_Collection,
            reactome_parser.PhysicalEntity_Collections)

        insert_collections(reactome_parser.Dna_Collections)
        insert_collections(reactome_parser.Rna_Collections)
        insert_collections(reactome_parser.SmallMolecule_Collections)
        insert_collections(reactome_parser.Protein_Collections)
        insert_collections(reactome_parser.Complex_Collections)
        insert_collections(reactome_parser.PhysicalEntity_Collections)

        insert_complex_parts(reactome_parser.Complexes)

        # NOW dump the ForbiddenIDs
        pickle.dump(ForbiddenIDs, open(Dumps.Forbidden_IDs, 'w'))

    if skip_import == 'M':
        get_all_meta_sets()

    # Meta insert/retrieval finished

    insert_reactions(DatabaseGraph.TemplateReaction, reactome_parser.TemplateReactions)
    insert_reactions(DatabaseGraph.Degradation, reactome_parser.Degradations)
    insert_reactions(DatabaseGraph.BiochemicalReaction, reactome_parser.BiochemicalReactions)

    # Reaction insert finished
    insert_catalysis(reactome_parser.Catalysises)
    insert_modulation(reactome_parser.Modulations)
    insert_pathways(reactome_parser.PathwaySteps, reactome_parser.Pathways)


if __name__ == "__main__":
    # insert_all()
    # run_diagnostics(bulbs_names_dict)
    # clear_all(bulbs_names_dict)
    pass
