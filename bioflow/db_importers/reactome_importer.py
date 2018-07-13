"""
Created on Jun 15, 2013
:author: andrei
"""
import pickle
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import Dumps, reactome_biopax_path, verbosity
from bioflow.internal_configs import Leg_ID_Filter
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph, on_alternative_graph
from bioflow.neo4j_db.db_io_routines import get_db_id
from bioflow.bio_db_parsers.reactomeParser import ReactomeParser
from bioflow.neo4j_db.db_io_routines import _bulb_specific_stable_get_all
from bioflow.neo4j_db.graph_content import neo4j_names_dict

log = get_logger(__name__)

# TODO: REFACTORING: put everything that uses memoization dictionary into a class
# TODO: with the new engine, memoization dictionary makes no more sense
memoization_dict = {}  # accelerated access pointer to the objects
ForbiddenIDs = []


def insert_cell_locations(cell_locations_dict):
    """
    Creates nodes corresponding to cell locations

    :param cell_locations_dict:
    """
    for Loc, displayName in cell_locations_dict.iteritems():
        if on_alternative_graph:
            memoization_dict[Loc] = DatabaseGraph.create('Location',
                                                         {'legacyId': Loc,
                                                          'displayName': displayName})
        else:
            memoization_dict[Loc] = DatabaseGraph.Location.create(ID=Loc, displayName=displayName)


def insert_minimal_annotations(annotated_node, annot_type_2_annot_list):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct,
    local DB ID and thus) needs to be inserted at the same time as the annotated object

    :param annotated_node:
    :param annot_type_2_annot_list:
    """
    if on_alternative_graph:
        DatabaseGraph.attach_all_node_annotations(annotated_node.id, annot_type_2_annot_list)
    else:
        for annotation_type, annotation_set in annot_type_2_annot_list.iteritems():
            if annotation_type != 'name' and annotation_set != '' and annotation_set != []:
                if not isinstance(annotation_type, list):
                    annotation_set = [annotation_set]
                for annotation in annotation_set:
                    annotation_node = DatabaseGraph.AnnotNode.create(
                        ptype=annotation_type, payload=annotation)
                    if verbosity > 1:
                        log.info('created annotation %s for %s, _id:%s',
                                 annotation_node, annotated_node, annotated_node.ID)
                    else:
                        log.debug('created annotation %s for %s, _id:%s',
                                 annotation_node, annotated_node, annotated_node.ID)
                    DatabaseGraph.is_annotated.create(
                        annotated_node,
                        annotation_node,
                        costum_from=annotated_node.ID,
                        costum_to='Annotation')


def insert_meta_objects(neo4j_graph_class, meta_id_2_property_dict):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a
     bulbs class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param meta_id_2_property_dict:
    """
    def bulbs_way():
        total_properties = len(meta_id_2_property_dict)

        config = neo4j_graph_class.client.config
        element_type = neo4j_graph_class.element_class.get_element_type(config)

        for i, (meta_name, property_dict) in enumerate(meta_id_2_property_dict.iteritems()):
            if i*20 % total_properties < 21:
                log.info('insert %s: %s out of %s', element_type, i, total_properties)

            primary = neo4j_graph_class.create(
                ID=meta_name,
                displayName=property_dict['displayName'],
                localization=property_dict['cellularLocation'],
                main_connex=False)

            if meta_name in Leg_ID_Filter:
                ForbiddenIDs.append(get_db_id(primary))

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
                            costum_from=primary.ID, costum_to=located_modification.ID,
                            source='Reactome_modification')

    def native_way():
        for i, (meta_name, property_dict) in enumerate(meta_id_2_property_dict.iteritems()):
            meta_properties = {'legacyId': meta_name,
                               'displayName': property_dict['displayName'],
                               'localization': property_dict['cellularLocation'],
                               'main_connex': False}
            primary = DatabaseGraph.create(neo4j_graph_class, meta_properties)

            if meta_name in Leg_ID_Filter:
                ForbiddenIDs.append(get_db_id(primary))

            memoization_dict[meta_name] = primary

            insert_minimal_annotations(primary, property_dict['references'])

            if 'cellularLocation' in property_dict.keys():
                secondary = memoization_dict[property_dict['cellularLocation']]
                DatabaseGraph.link(primary.id,
                                   secondary.id,
                                   'is_localized',
                                   {'custom_from': primary.properties['legacyId'],
                                    'custom_to': secondary.properties['legacyId']})

            if 'modification' in property_dict.keys():
                for modification in property_dict['modification']:
                    if 'location' in modification.keys() and 'modification' in modification.keys():

                        located_modification = DatabaseGraph.create('ModificationFeature',
                                                                    {'legacyId': modification['ID'],
                                                                     'type': 'post-translational_Mod',
                                                                     'location': modification['location'],
                                                                     'displayName': modification['modification']})

                        DatabaseGraph.link(primary.id,
                                           located_modification.id,
                                           'is_able_to_modify',
                                           {'custom_from': primary.properties['legacyId'],
                                            'custom_to': located_modification.properties['legacyId'],
                                            'source': 'Reactome_modification'})

    if on_alternative_graph:
        native_way()
    else:
        bulbs_way()


def insert_collections(collections_2_members):
    """
    Links a collection object reference to the members of the collection.

    :param collections_2_members:
    """
    for collection, collection_property_dict in collections_2_members.iteritems():
        # TODO: add the display name
        for member in collection_property_dict['collectionMembers']:
            if on_alternative_graph:
                collection_node = memoization_dict[collection]
                member_node = memoization_dict[member]
                DatabaseGraph.link(collection_node.id,
                                   member_node.id,
                                   'is_part_of_collection',
                                   {'custom_from': collection_node.properties['legacyId'],
                                    'custom_to': member_node.properties['legacyId'],
                                    'source': 'Reactome_collection'})
            else:
                DatabaseGraph.is_part_of_collection.create(
                    memoization_dict[collection],
                    memoization_dict[member],
                    costum_from=memoization_dict[collection].ID,
                    costum_to=memoization_dict[member].ID,
                    source='Reactome_collection')


def insert_complex_parts(complex_property_dict):
    """
    Links part of a complex to the complex

    :param complex_property_dict:
    """
    for key in complex_property_dict.keys():
        for part in complex_property_dict[key]['parts']:
            # TODO: remove redundant protection from Stoichiometry
            if 'Stoichiometry' not in part:
                if on_alternative_graph:
                    complex_node = memoization_dict[key]
                    part_node = memoization_dict[part]
                    DatabaseGraph.link(complex_node.id, part_node.id, 'is_part_of_complex',
                                       {'custom_from': complex_node.properties['legacyId'],
                                        'custom_to': part_node.properties['legacyId'],
                                        'source': 'Reactome_complex'})
                else:
                    DatabaseGraph.is_part_of_complex.create(
                        memoization_dict[key],
                        memoization_dict[part],
                        costum_from=memoization_dict[key].ID,
                        costum_to=memoization_dict[part].ID,
                        source='Reactome_complex')


def insert_reactions(neo4j_graph_class, property_source_dict):
    """
    Inserts a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs
    class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param property_source_dict:
    """
    for reaction, reaction_properties in property_source_dict.iteritems():
        if on_alternative_graph:
            memoization_dict[reaction] = DatabaseGraph.create(neo4j_graph_class,
                                                              {'legacyId': reaction,
                                                               'displayName': reaction_properties['displayName']})
        else:
            memoization_dict[reaction] = neo4j_graph_class.create(
                ID=reaction, displayName=reaction_properties['displayName'])

        insert_minimal_annotations(
            memoization_dict[reaction],
            reaction_properties['references'])

        for property_name, property_value_list in reaction_properties.iteritems():
            if property_name in ['left', 'right']:
                for elt in property_value_list:
                    if on_alternative_graph:
                        reaction_node = memoization_dict[reaction]
                        elt_node = memoization_dict[elt]
                        DatabaseGraph.link(reaction_node.id,
                                           elt_node.id,
                                           'is_reaction_participant',
                                           {'side': property_name,
                                            'custom_from': reaction_node.properties['legacyId'],
                                            'custom_to': reaction_node.properties['legacyId'],
                                            'source': 'Reactome_reaction'})
                    else:
                        DatabaseGraph.is_reaction_participant.create(
                            memoization_dict[reaction],
                            memoization_dict[elt],
                            side=property_name,
                            costum_from=memoization_dict[reaction].ID,
                            costum_to=memoization_dict[elt].ID,
                            source='Reactome_reaction')


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

                if on_alternative_graph:
                    memoization_dict[catalysis] = DatabaseGraph.link(
                        controller.id,
                        controlled.id,
                        'is_catalysant',
                        {'legacyId': catalysis,
                            'controlType': catalysis_properties['ControlType'],
                            'custom_from': controller.properties['legacyId'],
                            'custom_to': controlled.properties['legacyId'],
                            'source': 'Reactome_catalysis'})
                else:
                    memoization_dict[catalysis] = DatabaseGraph.is_catalysant.create(  #
                        controller,
                        controlled,
                        ID=catalysis,
                        controlType=catalysis_properties['ControlType'],
                        costum_from=controller.ID,
                        costum_to=controlled.ID,
                        source='Reactome_catalysis')

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
        if on_alternative_graph:
            memoization_dict[modulation] = DatabaseGraph.link(
                controller.id,
                controlled.id,
                'is_regulant',
                {'legacyID': modulation,
                 'controlType': modulation_property_dict['controlType'],
                 'custom_from': controller.properties['legacyId'],
                 'custom_to': controlled.properties['legacyId']
                 }
            )
        else:
            memoization_dict[modulation] = DatabaseGraph.is_regulant.create(       #
                controller,
                controlled,
                ID=modulation,
                controlType=modulation_property_dict['controlType'],
                costum_from=controller.ID,
                costum_to=controlled.ID,
                source='Reactome_modulation')


def insert_pathways(pathway_steps, pathways):
    """
    Inserts all the Pathways, linking and chaining subpathways
    Attention, it have to be imported at the same time as the reactions.

    :param pathway_steps:
    :param pathways:
    """
    for pathway_step in pathway_steps.keys():
        if on_alternative_graph:
            memoization_dict[pathway_step] = DatabaseGraph.create('PathwayStep',
                                                                  {'legacyId': pathway_step})
        else:
            memoization_dict[pathway_step] = DatabaseGraph.PathwayStep.create(ID=pathway_step)

    for pathway_step in pathways.keys():
        if on_alternative_graph:
            memoization_dict[pathway_step] = DatabaseGraph.create('Pathway',
                                                                  {'legacyId': pathway_step,
                                                                   'displayName': pathways[pathway_step]['displayName']})
        else:
            memoization_dict[pathway_step] = DatabaseGraph.Pathway.create(
                ID=pathway_step, displayName=pathways[pathway_step]['displayName'])

    for pathway_step in pathway_steps.keys():

        for component in pathway_steps[pathway_step]['components']:
            if on_alternative_graph:
                DatabaseGraph.link(memoization_dict[pathway_step].id,
                                   memoization_dict[component].id,
                                   'is_part_of_pathway',
                                   {'custom_from': pathway_step,
                                    'custom_to': component,
                                    'source': 'Reactome_pathway'})
            else:
                DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                        memoization_dict[component],
                                                        costum_from=pathway_step,
                                                        costum_to=component,
                                                        source='Reactome_pathway')

        for next_step in pathway_steps[pathway_step]['nextStep']:
            if on_alternative_graph:
                DatabaseGraph.link(memoization_dict[pathway_step].id,
                                   memoization_dict[next_step].id,
                                   'is_next_in_pathway',
                                   {'custom_from': pathway_step,
                                    'custom_to': next_step,
                                    'source': 'Reactome_pathway'})
            else:
                DatabaseGraph.is_next_in_pathway.create(memoization_dict[pathway_step],
                                                        memoization_dict[next_step],
                                                        costum_from=pathway_step,
                                                        costum_to=next_step,
                                                        source='Reactome_pathway')

    for pathway_step in pathways.keys():

        for second_pathway_step in pathways[pathway_step]['PathwayStep']:
            if on_alternative_graph:
                DatabaseGraph.link(memoization_dict[pathway_step].id,
                                   memoization_dict[second_pathway_step].id,
                                   'is_part_of_pathway',
                                   {'custom_from': pathway_step,
                                    'custom_to': second_pathway_step,
                                    'source': 'Reactome_pathway'})
            else:
                DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                        memoization_dict[second_pathway_step],
                                                        costum_from=pathway_step,
                                                        costum_to=second_pathway_step,
                                                        source='Reactome_pathway')

        for sub_pathway in pathways[pathway_step]['components']:
            if on_alternative_graph:
                DatabaseGraph.link(memoization_dict[pathway_step].id,
                                   memoization_dict[sub_pathway].id,
                                   'is_part_of_pathway',
                                   {'custom_from': pathway_step,
                                    'custom_to': sub_pathway,
                                    'source': 'Reactome_pathway'})
            else:
                DatabaseGraph.is_part_of_pathway.create(memoization_dict[pathway_step],
                                                        memoization_dict[sub_pathway],
                                                        costum_from=pathway_step,
                                                        costum_to=sub_pathway,
                                                        source='Reactome_pathway')


def get_one_meta_set(neo4j_class):
    """
    In case a MetaObject was already inserted, reloads it to the local dictionary for further
    annotation insertion

    :param neo4j_class:
    """
    if on_alternative_graph:
        for node in DatabaseGraph.get_all(neo4j_class):
            memoization_dict[node.id] = node
    else:
        for bulbs_object in _bulb_specific_stable_get_all(neo4j_class):
            if bulbs_object is not None:
                memoization_dict[bulbs_object.ID] = bulbs_object


def get_all_meta_sets():
    """
    In case the MetaObjects were already inserted, reloads them all to the local dictionary for
    further annotation insertion
    """

    if on_alternative_graph:
        list_of_bulbs_classes = [
            neo4j_names_dict['DNA'],
            neo4j_names_dict['DNA Collection'],
            neo4j_names_dict['RNA'],
            neo4j_names_dict['RNA Collection'],
            neo4j_names_dict['Small Molecule'],
            neo4j_names_dict['Small Molecule Collection'],
            neo4j_names_dict['Protein'],
            neo4j_names_dict['Protein Collection'],
            neo4j_names_dict['Complex'],
            neo4j_names_dict['Complex Collection'],
            neo4j_names_dict['Physical Entity'],
            neo4j_names_dict['Physical Entity Collection']
        ]
    else:
        list_of_bulbs_classes = [
            neo4j_names_dict['DNA'][0],
            neo4j_names_dict['DNA Collection'][0],
            neo4j_names_dict['RNA'][0],
            neo4j_names_dict['RNA Collection'][0],
            neo4j_names_dict['Small Molecule'][0],
            neo4j_names_dict['Small Molecule Collection'][0],
            neo4j_names_dict['Protein'][0],
            neo4j_names_dict['Protein Collection'][0],
            neo4j_names_dict['Complex'][0],
            neo4j_names_dict['Complex Collection'][0],
            neo4j_names_dict['Physical Entity'][0],
            neo4j_names_dict['Physical Entity Collection'][0]
        ]

    for node_class in list_of_bulbs_classes:
        get_one_meta_set(node_class)


def insert_reactome(skip_import='N'):
    """
    Performs the massive import of the Reactome database into the local neo4j database.

    :param skip_import:     * N => will skip nothing and implement the import once and for all.
                     * M => skips meta import, recovers the metas and resumes from the Reactions
                     import.
    """
    def get_db_class_handle(shortname):
        if on_alternative_graph:
            return neo4j_names_dict[shortname]
        else:
            return neo4j_names_dict[shortname][0]

    reactome_parser = ReactomeParser(reactome_biopax_path)
    reactome_parser.parse_all()

    if skip_import == 'N':

        insert_cell_locations(reactome_parser.CellularLocations)

        insert_meta_objects(get_db_class_handle('DNA'), reactome_parser.Dnas)
        insert_meta_objects(get_db_class_handle('DNA Collection'), reactome_parser.Dna_Collections)
        insert_meta_objects(get_db_class_handle('RNA'), reactome_parser.Rnas)
        insert_meta_objects(get_db_class_handle('RNA Collection'), reactome_parser.Rna_Collections)
        insert_meta_objects(get_db_class_handle('Small Molecule'), reactome_parser.SmallMolecules)
        insert_meta_objects(
            get_db_class_handle('Small Molecule Collection'),
            reactome_parser.SmallMolecule_Collections)
        insert_meta_objects(get_db_class_handle('Protein'), reactome_parser.Proteins)
        insert_meta_objects(get_db_class_handle('Protein Collection'), reactome_parser.Protein_Collections)
        insert_meta_objects(get_db_class_handle('Complex'), reactome_parser.Complexes)
        insert_meta_objects(get_db_class_handle('Complex Collection'), reactome_parser.Complex_Collections)
        insert_meta_objects(get_db_class_handle('Physical Entity'), reactome_parser.PhysicalEntities)
        insert_meta_objects(
            get_db_class_handle('Physical Entity Collection'),
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


    # TODO: add insertion log to prevent the impression of hanging
    if skip_import == 'M':
        get_all_meta_sets()

    # Meta insert/retrieval finished

    insert_reactions(get_db_class_handle('Template Reaction'), reactome_parser.TemplateReactions)
    insert_reactions(get_db_class_handle('Degradation'), reactome_parser.Degradations)
    insert_reactions(get_db_class_handle('BiochemicalReaction'), reactome_parser.BiochemicalReactions)

    # Reaction insert finished
    insert_catalysis(reactome_parser.Catalysises)
    insert_modulation(reactome_parser.Modulations)
    insert_pathways(reactome_parser.PathwaySteps, reactome_parser.Pathways)


if __name__ == "__main__":
    # insert_all()
    # run_diagnostics(neo4j_names_dict)
    # clear_all(neo4j_names_dict)
    pass
