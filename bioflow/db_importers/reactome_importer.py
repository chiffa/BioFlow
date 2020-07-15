"""
Created on Jun 15, 2013
:author: andrei
"""
import pickle
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import Dumps, reactome_biopax_path
from bioflow.internal_configs import Leg_ID_Filter, neo4j_names_dict
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.neo4j_db.db_io_routines import get_db_id
from bioflow.bio_db_parsers.reactomeParser import ReactomeParser

log = get_logger(__name__)

# TODO: REFACTORING: put everything that uses memoization dictionary into a class

memoization_dict = {}  # accelerated access pointer to the objects
ForbiddenIDs = []


def insert_cell_locations(cell_locations_dict):
    """
    Creates nodes corresponding to cell locations

    :param cell_locations_dict:
    """
    for Loc, displayName in cell_locations_dict.items():
        memoization_dict[Loc] = DatabaseGraph.create('Location',
                                                     {'legacyId': Loc,
                                                      'displayName': displayName})


def insert_minimal_annotations(annotated_node, annot_type_2_annot_list):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct,
    local DB ID and thus) needs to be inserted at the same time as the annotated object

    :param annotated_node:
    :param annot_type_2_annot_list:
    """
    DatabaseGraph.attach_all_node_annotations(annotated_node.id, annot_type_2_annot_list)


def insert_meta_objects(neo4j_graph_class, meta_id_2_property_dict):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a
     bulbs class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param meta_id_2_property_dict:
    """
    size = len(list(meta_id_2_property_dict.keys()))
    log.info('Starting inserting %s with %s elements', neo4j_graph_class, size)
    breakpoints = 300
    for i, (meta_name, property_dict) in enumerate(meta_id_2_property_dict.items()):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i)/float(size)*100.0))

        meta_properties = {'legacyId': meta_name,
                           'displayName': property_dict['displayName'],
                           'localization': property_dict['cellularLocation'],
                           'main_connex': False}
        primary = DatabaseGraph.create(neo4j_graph_class, meta_properties)
        print(primary)
        print(dir(primary))
        print(primary._properties, '\n',
              primary.values(), '\n',
              primary.labels, '\n',
              primary.keys(), '\n',
              primary.items())

        if meta_name in Leg_ID_Filter:
            ForbiddenIDs.append(get_db_id(primary))

        memoization_dict[meta_name] = primary

        insert_minimal_annotations(primary, property_dict['references'])

        if 'cellularLocation' in list(property_dict.keys()):
            secondary = memoization_dict[property_dict['cellularLocation']]
            DatabaseGraph.link(primary.id,
                               secondary.id,
                               'is_localized',
                               {'custom_from': primary._properties['legacyId'],
                                'custom_to': secondary._properties['legacyId']})

        if 'modification' in list(property_dict.keys()):
            for modification in property_dict['modification']:
                if 'location' in list(modification.keys()) and 'modification' in list(modification.keys()):

                    located_modification = DatabaseGraph.create('ModificationFeature',
                                                                {'legacyId': modification['ID'],
                                                                 'type': 'post-translational_Mod',
                                                                 'location': modification['location'],
                                                                 'displayName': modification['modification']})

                    DatabaseGraph.link(primary.id,
                                       located_modification.id,
                                       'is_able_to_modify',
                                       {'custom_from': primary._properties['legacyId'],
                                        'custom_to': located_modification._properties['legacyId'],
                                        'source': 'Reactome_modification'})


def insert_collections(collections_2_members):
    """
    Links a collection object reference to the members of the collection.

    :param collections_2_members:
    """

    breakpoints = 300
    size = len(list(collections_2_members.keys()))

    for i, (collection, collection_property_dict) in enumerate(collections_2_members.items()):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i)/float(size)*100))


        for member in collection_property_dict['collectionMembers']:
            collection_node = memoization_dict[collection]
            member_node = memoization_dict[member]
            DatabaseGraph.link(collection_node.id,
                               member_node.id,
                               'is_part_of_collection',
                               {'custom_from': collection_node.properties['legacyId'],
                                'custom_to': member_node.properties['legacyId'],
                                'source': 'Reactome_collection'})


def insert_complex_parts(complex_property_dict):
    """
    Links part of a complex to the complex

    :param complex_property_dict:
    """
    breakpoint = 300
    size = len(list(complex_property_dict.keys()))

    for i, key in enumerate(complex_property_dict.keys()):

        if i % breakpoint == 0:
            log.info('\t %.2f %%' % (float(i)/float(size)*100.0))

        for part in complex_property_dict[key]['parts']:
            # TODO: remove redundant protection from Stoichiometry
            if 'Stoichiometry' not in part:
                complex_node = memoization_dict[key]
                part_node = memoization_dict[part]
                DatabaseGraph.link(complex_node.id, part_node.id, 'is_part_of_complex',
                                   {'custom_from': complex_node.properties['legacyId'],
                                    'custom_to': part_node.properties['legacyId'],
                                    'source': 'Reactome_complex'})


def insert_reactions(neo4j_graph_class, property_source_dict):
    """
    Inserts a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs
    class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param property_source_dict:
    """
    for reaction, reaction_properties in property_source_dict.items():
        memoization_dict[reaction] = DatabaseGraph.create(neo4j_graph_class,
                                                          {'legacyId': reaction,
                                                           'displayName': reaction_properties['displayName']})

        insert_minimal_annotations(
            memoization_dict[reaction],
            reaction_properties['references'])

        for property_name, property_value_list in reaction_properties.items():
            if property_name in ['left', 'right']:
                for elt in property_value_list:
                    reaction_node = memoization_dict[reaction]
                    elt_node = memoization_dict[elt]
                    DatabaseGraph.link(reaction_node.id,
                                       elt_node.id,
                                       'is_reaction_participant',
                                       {'side': property_name,
                                        'custom_from': reaction_node.properties['legacyId'],
                                        'custom_to': reaction_node.properties['legacyId'],
                                        'source': 'Reactome_reaction'})


# TODO: catalysis and modulation are identical in the marked lines
def insert_catalysis(catalysises_dict):
    """
    Inserts all the catalysis links from one meta-element to an another

    :param catalysises_dict:
    """
    for catalysis, catalysis_properties in catalysises_dict.items():

        if 'controller' in list(catalysis_properties.keys()) \
                and 'controlled' in list(catalysis_properties.keys()):

            if catalysis_properties['controlled'] in list(memoization_dict.keys()) \
                    and catalysis_properties['controller'] in list(memoization_dict.keys()):

                # TODO: this can be moved to the parsing
                if 'ControlType' not in list(catalysises_dict[catalysis].keys()):
                    catalysis_properties['ControlType'] = 'UNKNOWN'

                controller = memoization_dict[catalysis_properties['controller']]  #
                controlled = memoization_dict[catalysis_properties['controlled']]  #

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
                log.debug("Catalysis targets not memoized: %s : %s, %s, %s", catalysis,
                          catalysises_dict[catalysis],
                          catalysises_dict[catalysis]['controlled'] in list(memoization_dict.keys()),
                          catalysises_dict[catalysis]['controller'] in list(memoization_dict.keys()))
        else:
            log.debug("Catalysis without control/controlled %s : %s, %s, %s,",
                      catalysis, catalysises_dict[catalysis],
                      'controller' in list(catalysises_dict[catalysis].keys()),
                      'controlled' in list(catalysises_dict[catalysis].keys()))


def insert_modulation(modulations_dict):
    """
    Inserts all the Modulation links from one meta-element to an another

    :param modulations_dict:
    """
    for modulation, modulation_property_dict in modulations_dict.items():
        controller = memoization_dict[modulation_property_dict['controller']]  #
        controlled = memoization_dict[modulation_property_dict['controlled']]  #

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


def insert_pathways(pathway_steps, pathways):
    """
    Inserts all the Pathways, linking and chaining subpathways
    Attention, it have to be imported at the same time as the reactions.

    :param pathway_steps:
    :param pathways:
    """
    breakpoints = 300
    ps_len = len(list(pathway_steps.keys()))
    p_len = len(list(pathways.keys()))

    log.info('Inserting Pathway steps with %s elements', len(list(pathway_steps.keys())))

    for i, pathway_step in enumerate(pathway_steps.keys()):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i)/float(ps_len)*100))

        memoization_dict[pathway_step] = DatabaseGraph.create('PathwayStep',
                                                              {'legacyId': pathway_step})

    log.info('Inserting Pathways with %s elements', len(list(pathways.keys())))

    for pathway_step in list(pathways.keys()):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i)/float(p_len)*100))


        memoization_dict[pathway_step] = DatabaseGraph.create('Pathway',
                                                              {'legacyId': pathway_step,
                                                               'displayName': pathways[pathway_step]['displayName']})

    for i, pathway_step in enumerate(pathway_steps.keys()):

        for component in pathway_steps[pathway_step]['components']:
            DatabaseGraph.link(memoization_dict[pathway_step].id,
                               memoization_dict[component].id,
                               'is_part_of_pathway',
                               {'custom_from': pathway_step,
                                'custom_to': component,
                                'source': 'Reactome_pathway'})

        for next_step in pathway_steps[pathway_step]['nextStep']:
            DatabaseGraph.link(memoization_dict[pathway_step].id,
                               memoization_dict[next_step].id,
                               'is_next_in_pathway',
                               {'custom_from': pathway_step,
                                'custom_to': next_step,
                                'source': 'Reactome_pathway'})

    for pathway_step in list(pathways.keys()):

        for second_pathway_step in pathways[pathway_step]['PathwayStep']:
            DatabaseGraph.link(memoization_dict[pathway_step].id,
                               memoization_dict[second_pathway_step].id,
                               'is_part_of_pathway',
                               {'custom_from': pathway_step,
                                'custom_to': second_pathway_step,
                                'source': 'Reactome_pathway'})

        for sub_pathway in pathways[pathway_step]['components']:
            DatabaseGraph.link(memoization_dict[pathway_step].id,
                               memoization_dict[sub_pathway].id,
                               'is_part_of_pathway',
                               {'custom_from': pathway_step,
                                'custom_to': sub_pathway,
                                'source': 'Reactome_pathway'})


def get_one_meta_set(neo4j_class):
    """
    In case a MetaObject was already inserted, reloads it to the local dictionary for further
    annotation insertion

    :param neo4j_class:
    """
    for node in DatabaseGraph.get_all(neo4j_class):
        memoization_dict[node.properties['legacyId']] = node


def get_all_meta_sets():
    """
    In case the MetaObjects were already inserted, reloads them all to the local dictionary for
    further annotation insertion
    """
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
        return neo4j_names_dict[shortname]

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

        log.info('Inserting DNA Collections with %s elements' % len(reactome_parser.Dna_Collections))
        insert_collections(reactome_parser.Dna_Collections)
        log.info('Inserting RNA Collections with %s elements' % len(reactome_parser.Rna_Collections))
        insert_collections(reactome_parser.Rna_Collections)
        log.info('Inserting Small Molecule Collections with %s elements' % len(reactome_parser.SmallMolecule_Collections))
        insert_collections(reactome_parser.SmallMolecule_Collections)
        log.info('Inserting Protein Collections with %s elements' % len(reactome_parser.Protein_Collections))
        insert_collections(reactome_parser.Protein_Collections)
        log.info('Inserting Complex Collections with %s elements' % len(reactome_parser.Complex_Collections))
        insert_collections(reactome_parser.Complex_Collections)
        log.info('Inserting Physical Entity Collections with %s elements' % len(reactome_parser.PhysicalEntity_Collections))
        insert_collections(reactome_parser.PhysicalEntity_Collections)

        log.info('Inserting Complexes with %s elements' % (len(reactome_parser.Complexes)))
        insert_complex_parts(reactome_parser.Complexes)

        # NOW dump the ForbiddenIDs
        pickle.dump(ForbiddenIDs, open(Dumps.Forbidden_IDs, 'wt'))


    if skip_import == 'M':
        get_all_meta_sets()

    # print memoization_dict.keys()

    # Meta insert/retrieval finished
    log.info('Inserting Template Reactions with %s elements' % len(reactome_parser.Dna_Collections))
    insert_reactions(get_db_class_handle('TemplateReaction'), reactome_parser.TemplateReactions)
    log.info('Inserting Degradations with %s elements' % len(reactome_parser.Dna_Collections))
    insert_reactions(get_db_class_handle('Degradation'), reactome_parser.Degradations)
    log.info('Inserting Biochemical Reactions with %s elements' % len(reactome_parser.Dna_Collections))
    insert_reactions(get_db_class_handle('BiochemicalReaction'), reactome_parser.BiochemicalReactions)

    # Reaction insert finished
    log.info('Inserting Catalyses with %s elements' % len(reactome_parser.Catalysises))
    insert_catalysis(reactome_parser.Catalysises)
    log.info('Inserting Modulations with %s elements' % len(reactome_parser.Modulations))
    insert_modulation(reactome_parser.Modulations)
    insert_pathways(reactome_parser.PathwaySteps, reactome_parser.Pathways)
