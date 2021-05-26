"""
Created on Jun 15, 2013
:author: andrei
"""
import pickle
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import Dumps, reactome_biopax_path
from bioflow.configs.main_configs import reactome_forbidden_nodes
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.bio_db_parsers.reactomeParser import ReactomeParser

log = get_logger(__name__)

memoization_dict = {}  # accelerated access pointer to the objects
ForbiddenIDs = []


def insert_cell_locations(cell_locations_dict):
    """
    Creates nodes corresponding to cell locations

    :param cell_locations_dict:
    """
    for Loc, displayName in cell_locations_dict.items():
        memoization_dict[Loc] = DatabaseGraph.create('Location',
                                                     {'legacyID': Loc,
                                                      'displayName': displayName,
                                                      'parse_type': 'annotation',
                                                      'source': 'Reactome'})


def insert_minimal_annotations(annotated_node, annot_type_2_annot_list, source):
    """
    Inserts a minimal annotation provided the annotated_node Node (it requires the direct,
    local DB ID and thus) needs to be inserted at the same time as the annotated object

    :param annotated_node:
    :param annot_type_2_annot_list:
    """
    DatabaseGraph.attach_all_node_annotations(annotated_node.id,
                                              annot_type_2_annot_list,
                                              source=source)


def insert_reactome_class(neo4j_graph_class, reactome_obj_id_2_property_dict, parse_type):
    """
    Inserst a Meta-Object (I.e. any physical entity or collection thereof) as a member of a
     bulbs class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param reactome_obj_id_2_property_dict:
    """
    size = len(list(reactome_obj_id_2_property_dict.keys()))
    log.info('Starting inserting %s with %s elements', neo4j_graph_class, size)
    breakpoints = 300

    for i, (reactome_id, property_dict) in enumerate(reactome_obj_id_2_property_dict.items()):

        if i % breakpoints == 0:
            # TODO: [progress bar]
            log.info('\t %.2f %%' % (float(i) / float(size) * 100.0))

        reactome_obj_properties = {'legacyID': reactome_id,
                                   'displayName': property_dict['displayName'],
                                   'localization':
                                       memoization_dict[property_dict['cellularLocation']]['displayName'],
                                   'source': 'Reactome',
                                   'parse_type': parse_type,
                                   'main_connex': False}

        primary = DatabaseGraph.create(neo4j_graph_class, reactome_obj_properties)

        log.debug(primary)
        log.debug(dir(primary))
        log.debug("%s\n%s\n%s\n%s\n%s" %
                  (primary._properties,
                  primary.values(),
                  primary.labels,
                  primary.keys(),
                  primary.items()))

        # print(primary)
        # print(dir(primary))
        # print(primary._properties, '\n',
        #       primary.values(), '\n',
        #       primary.labels, '\n',
        #       primary.keys(), '\n',
        #       primary.items())

        if reactome_id in reactome_forbidden_nodes:
            ForbiddenIDs.append(primary.id)

        memoization_dict[reactome_id] = primary

        insert_minimal_annotations(primary,
                                   property_dict['references'],
                                   source='Reactome')

        if 'cellularLocation' in list(property_dict.keys()):
            secondary = memoization_dict[property_dict['cellularLocation']]
            DatabaseGraph.link(primary.id,
                               secondary.id,
                               'is_localized',
                               {'source': 'Reactome',
                                'parse_type': 'annotates'})

        if 'modification' in list(property_dict.keys()):

            for modification in property_dict['modification']:

                if 'location' in list(modification.keys()) and 'modification' in list(modification.keys()):

                    located_modification = DatabaseGraph.create('ModificationFeature',
                                                                {'legacyID': modification['ID'],
                                                                 'type': 'post-translational_Mod',
                                                                 'location': modification['location'],
                                                                 'displayName': modification['modification'],
                                                                 'source': 'Reactome',
                                                                 'parse_type': 'physical_entity'})

                    DatabaseGraph.link(primary.id,
                                       located_modification.id,
                                       'is_able_to_modify',
                                       {'source_note': 'Reactome_modification',
                                        'source': 'Reactome',
                                        'parse_type': 'refines'}
                                       )


def insert_collections(collections_2_members):
    """
    Links a collection object reference to the members of the collection.

    :param collections_2_members:
    """

    breakpoints = 300
    size = len(list(collections_2_members.keys()))

    for i, (collection, collection_property_dict) in enumerate(collections_2_members.items()):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i) / float(size) * 100))


        for member in collection_property_dict['collectionMembers']:
            collection_node = memoization_dict[collection]
            member_node = memoization_dict[member]
            DatabaseGraph.link(collection_node.id,
                               member_node.id,
                               'is_part_of_collection',
                               {'source': 'Reactome',
                                'parse_type': 'refines',
                                'source_note': 'Reactome_collection'
                                })


def insert_complex_parts(complex_property_dict):
    """
    Links part of a complex to the complex

    :param complex_property_dict:
    """
    breakpoint = 300
    size = len(list(complex_property_dict.keys()))

    for i, key in enumerate(complex_property_dict.keys()):

        if i % breakpoint == 0:
            log.info('\t %.2f %%' % (float(i) / float(size) * 100.0))

        for part in complex_property_dict[key]['parts']:
            if 'Stoichiometry' not in part:
                complex_node = memoization_dict[key]
                part_node = memoization_dict[part]
                DatabaseGraph.link(complex_node.id, part_node.id, 'is_part_of_complex',
                                   {'source': 'Reactome',
                                    'parse_type': 'physical_entity_molecular_interaction',
                                    'source_note': 'Reactome_complex'
                                    })


def insert_reactions(neo4j_graph_class, property_source_dict):
    """
    Inserts a Reaction-Object (I.e. any reaction or type of reactions) as a member of a bulbs
    class and pumping the bioflow information from the property bioflow

    :param neo4j_graph_class:
    :param property_source_dict:
    """
    for reaction, reaction_properties in property_source_dict.items():
        memoization_dict[reaction] = DatabaseGraph.create(neo4j_graph_class,
                                                          {'legacyID': reaction,
                                                           'displayName': reaction_properties['displayName'],
                                                           'source': 'Reactome',
                                                           'parse_type': 'physical_entity'})

        insert_minimal_annotations(
            memoization_dict[reaction],
            reaction_properties['references'],
            source='Reactome')

        for property_name, property_value_list in reaction_properties.items():
            if property_name in ['left', 'right']:
                for elt in property_value_list:
                    reaction_node = memoization_dict[reaction]
                    elt_node = memoization_dict[elt]
                    DatabaseGraph.link(reaction_node.id,
                                       elt_node.id,
                                       'is_reaction_participant',
                                       {'side': property_name,
                                        'source_note': 'Reactome_reaction',
                                        'source': 'Reactome',
                                        'parse_type': 'physical_entity_molecular_interaction'})


# TODO: [data organization] catalysis need to be inserted as nodes and then cross-linked,
#  for better compatibility with the Pathway searcj
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

                if 'ControlType' not in list(catalysises_dict[catalysis].keys()):
                    catalysis_properties['ControlType'] = 'UNKNOWN'

                controller = memoization_dict[catalysis_properties['controller']]  #
                controlled = memoization_dict[catalysis_properties['controlled']]  #

                memoization_dict[catalysis] = DatabaseGraph.link(
                    controller.id,
                    controlled.id,
                    'is_catalysant',
                    {'legacyID': catalysis,
                        'controlType': catalysis_properties['ControlType'],
                        'source': 'Reactome',
                        'source_note': 'Reactome_catalysis',
                        'parse_type': 'physical_entity_molecular_interaction'
                     })

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
             'source': 'Reactome',
             'source_note': 'Reactome_modulation',
             'parse_type': 'physical_entity_molecular_interaction'
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
            log.info('\t %.2f %%' % (float(i) / float(ps_len) * 100))

        memoization_dict[pathway_step] = DatabaseGraph.create(
            'PathwayStep',
            {'legacyID': pathway_step,
             'displayName': pathway_steps[pathway_step].get('displayName', pathway_step),
             'source': 'Reactome',
             'parse_type': 'annotation'})

    log.info('Inserting Pathways with %s elements', len(list(pathways.keys())))


    # TODO: [reactome pathways sanity] links are inverted
    for i, pathway in enumerate(list(pathways.keys())):

        if i % breakpoints == 0:
            log.info('\t %.2f %%' % (float(i) / float(p_len) * 100))

        # print(pathways[pathway])

        memoization_dict[pathway] = DatabaseGraph.create(
            'Pathway',
            {'legacyID': pathway,
             'displayName': pathways[pathway].get('displayName', pathway),
             'source': 'Reactome',
             'parse_type': 'annotation'})

    for i, pathway_step in enumerate(pathway_steps.keys()):

        for component in pathway_steps[pathway_step]['components']:

            parse_type = 'annotation_relationship'
            if memoization_dict[component]['parse_type'] == 'physical_entity':
                parse_type = 'annotates'

            DatabaseGraph.link(memoization_dict[component].id,
                               memoization_dict[pathway_step].id,
                               'is_part_of_pathway',
                               {'source_note': 'Reactome_pathway',
                                'source': 'Reactome',
                                'parse_type': parse_type
                                })

        for next_step in pathway_steps[pathway_step]['nextStep']:
            # only links pathway steps
            DatabaseGraph.link(memoization_dict[pathway_step].id,
                               memoization_dict[next_step].id,
                               'is_next_in_pathway',
                               {'source_note': 'Reactome_pathway',
                                'source': 'Reactome',
                                'parse_type': 'annotation_relationship'
                                })

    for pathway in list(pathways.keys()):

        for second_pathway in pathways[pathway]['PathwayStep']:
            # only links to pathway steps
            DatabaseGraph.link(memoization_dict[second_pathway].id,
                               memoization_dict[pathway].id,
                               'is_part_of_pathway',
                               {'source_note': 'Reactome_pathway',
                                'source': 'Reactome',
                                'parse_type': 'annotation_relationship'
                                })

        for sub_pathway in pathways[pathway]['components']:

            parse_type = 'annotation_relationship'
            if memoization_dict[sub_pathway]['parse_type'] == 'physical_entity':
                parse_type = 'annotates'

            DatabaseGraph.link(memoization_dict[sub_pathway].id,
                               memoization_dict[pathway].id,
                               'is_part_of_pathway',
                               {'source_note': 'Reactome_pathway',
                                'source': 'Reactome',
                                'parse_type': parse_type
                                })


def re_memoize_reactome_nodes():
    """
    In case the Reactome meta objects were already inserted, reloads them all to the local
    dictionary for further annotation insertion
    """

    reactome_nodes = DatabaseGraph.find(filter_dict={'source': 'Reactome'})
    memoization_dict.update({node['legacyID']: node for node in reactome_nodes})


def insert_reactome(skip_import='N'):
    """
    Performs the massive import of the Reactome database into the local neo4j database.

    :param skip_import:     * N => will skip nothing and implement the import once and for all.
                     * M => skips meta import, recovers the metas and resumes from the Reactions
                     import.
    """
    reactome_parser = ReactomeParser(reactome_biopax_path)
    reactome_parser.parse_all()

    if skip_import == 'N':

        insert_cell_locations(reactome_parser.CellularLocations)

        insert_reactome_class('DNA',
                              reactome_parser.Dnas, 'physical_entity')
        insert_reactome_class("DNA_Collection",
                              reactome_parser.Dna_Collections, 'physical_entity')
        insert_reactome_class("RNA",
                              reactome_parser.Rnas, 'physical_entity')
        insert_reactome_class("RNA_Collection",
                              reactome_parser.Rna_Collections, 'physical_entity')
        insert_reactome_class("SmallMolecule",
                              reactome_parser.SmallMolecules, 'physical_entity')
        insert_reactome_class("SmallMolecule_Collection",
                              reactome_parser.SmallMolecule_Collections, 'physical_entity')
        insert_reactome_class("Protein",
                              reactome_parser.Proteins, 'physical_entity')
        insert_reactome_class("Protein_Collection",
                              reactome_parser.Protein_Collections, 'physical_entity')
        insert_reactome_class("Complex",
                              reactome_parser.Complexes, 'physical_entity')
        insert_reactome_class("Complex_Collection",
                              reactome_parser.Complex_Collections, 'physical_entity')
        insert_reactome_class("PhysicalEntity",
                              reactome_parser.PhysicalEntities, 'physical_entity')
        insert_reactome_class("PhysicalEntity_Collection",
                              reactome_parser.PhysicalEntity_Collections, 'physical_entity')

        log.info('Inserting DNA Collections with %s elements'
                 % len(reactome_parser.Dna_Collections))
        insert_collections(reactome_parser.Dna_Collections)
        log.info('Inserting RNA Collections with %s elements'
                 % len(reactome_parser.Rna_Collections))
        insert_collections(reactome_parser.Rna_Collections)
        log.info('Inserting Small Molecule Collections with %s elements'
                 % len(reactome_parser.SmallMolecule_Collections))
        insert_collections(reactome_parser.SmallMolecule_Collections)
        log.info('Inserting Protein Collections with %s elements'
                 % len(reactome_parser.Protein_Collections))
        insert_collections(reactome_parser.Protein_Collections)
        log.info('Inserting Complex Collections with %s elements'
                 % len(reactome_parser.Complex_Collections))
        insert_collections(reactome_parser.Complex_Collections)
        log.info('Inserting Physical Entity Collections with %s elements'
                 % len(reactome_parser.PhysicalEntity_Collections))
        insert_collections(reactome_parser.PhysicalEntity_Collections)

        log.info('Inserting Complexes with %s elements'
                 % (len(reactome_parser.Complexes)))
        insert_complex_parts(reactome_parser.Complexes)


    if skip_import == 'M':
        re_memoize_reactome_nodes()

    # print memoization_dict.keys()

    # Meta insert/retrieval finished
    log.info('Inserting Template Reactions with %s elements'
             % len(reactome_parser.TemplateReactions))
    insert_reactions("Template_Reaction",
                     reactome_parser.TemplateReactions)
    log.info('Inserting Degradations with %s elements'
             % len(reactome_parser.Degradations))
    insert_reactions("Degradation",
                     reactome_parser.Degradations)
    log.info('Inserting Biochemical Reactions with %s elements'
             % len(reactome_parser.BiochemicalReactions))
    insert_reactions("BiochemicalReaction",
                     reactome_parser.BiochemicalReactions)

    # Reaction insert finished
    log.info('Inserting Catalyses with %s elements'
             % len(reactome_parser.Catalysises))
    insert_catalysis(reactome_parser.Catalysises)
    # log.info('Inserting Modulations with %s elements' # ceased to exist
    #          % len(reactome_parser.Modulations))
    # insert_modulation(reactome_parser.Modulations)
    insert_pathways(reactome_parser.PathwaySteps,
                    reactome_parser.Pathways)

    # Q: There is no linking from the actual entities to pathways in Reactome?
    #   There are, but through pathway steps.
