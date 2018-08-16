"""
Module responsible for the declaration of all the database access routines that are used by other
modules. In case a different back-end than neo4j is used, all methods in this cluster have to be
re-implemented
"""
import os
import pickle
from collections import defaultdict
from csv import reader, writer
from pprint import PrettyPrinter, pprint
from bioflow.main_configs import forbidden_neo4j_ids, Dumps, analysis_protein_ids_csv, \
    background_protein_ids_csv, verbosity
from bioflow.internal_configs import edge_type_filters, Leg_ID_Filter, annotation_nodes_ptypes
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.neo4j_db.graph_content import neo4j_names_dict, full_list, forbidden_verification_list
from bioflow.utils.log_behavior import get_logger
from bioflow.utils.io_routines import write_to_csv, memoize, time_exection

log = get_logger(__name__)


def get_db_id(neo4j_node):
    """
    gets bulbs_id for a node, since it's a hidden attribute

    :param neo4j_node:
    :return:
    """
    return neo4j_node.id


# TODO: REFACTOR: offload annotations to ElasticSearch engine
def get_attached_annotations(neo4j_node_id):
    """
    Recovers ids of annotation nodes attached the the node with a given bulbs id

    :param neo4j_node_id:
    :return:
    """
    list_of_annotations = []
    annotation_node_generator = DatabaseGraph.get_linked(neo4j_node_id, link_type="is_annotated")
    if not annotation_node_generator:
        log.debug("node %s has not annotations attached to it", neo4j_node_id)
        return []
    else:
        for rel_node in annotation_node_generator:
            list_of_annotations.append(get_db_id(rel_node))
        return list_of_annotations


def node_generator_2_db_ids(node_generator):
    """
    Get bulb ids of nodes in the generator

    :param node_generator: iterator over annot_nodes
    :return: the bubls ids of nodes in generator
    """
    if not node_generator:
        return []

    node_set = []
    for node in node_generator:
        node_set.append(get_db_id(node))
        # print 'forbidding ', node.ID, node.displayName

    return node_set


def look_up_annotation_node(p_load, p_type=''):
    """
    Looks up nodes accessible via the annotation nodes with a given annotation and given
    annotation  type.
    The lookup strict match, but case-insensitive.

    :param p_load: payload
    :param p_type: payload type. servers to restrict search to a specific ID subset
    :return: (node type, node's displayName, node's db_ID, node's legacy ID)
    :raise Exception: "p_type unsupported", in case a p_type is not on the supported list
    specified  in the neo4j_db.neo4j_typeDec
    """

    payload = p_load.upper()  # make payload case-agnostic

    nodes = DatabaseGraph.get_from_annotation_tag(payload, p_type)
    retlist = []
    for node in nodes:
        node_bulbs_id = get_db_id(node)
        node_legacy_id = node.properties['legacyId']
        node_type = list(node.labels)[0]
        node_display_name = node.properties['displayName']
        retlist.append((node_type, node_display_name, node_bulbs_id, node_legacy_id))
    return retlist


def look_up_annotation_set(p_load_list, p_type=''):
    """
    Looks up an set of annotations in the database and finds the Ids of nodes containing SWISSPROT
    proteins linked to by annotations

    :param p_load_list: list of payloads
    :param p_type: expected type of payloads
    :return:
    """
    def db_id_mapping_helper(mapped_db_id_list):
        if mapped_db_id_list:
            return mapped_db_id_list[0][2]
        else:
            return ''

    def transform_annotation_node(neo4j_native_nodes):
        retlist = []
        for node in neo4j_native_nodes:
            node_bulbs_id = get_db_id(node)
            node_legacy_id = node.properties['legacyId']
            node_type = list(node.labels)[0]
            node_display_name = node.properties['displayName']
            payload = (node_type, node_display_name, node_bulbs_id, node_legacy_id)
            if node_type == 'UNIPROT':
                retlist.insert(0, payload)
            else:
                retlist.append(payload)
        return retlist

    load_2_name_list = [(p_load, transform_annotation_node(p_nodes)) for (p_load, p_nodes) in
                        zip(p_load_list, DatabaseGraph.batch_retrieve_from_annotation_tags(p_load_list, p_type))]

    db_id_list = [db_id_mapping_helper(value) for key, value in load_2_name_list]

    not_found_list = [key for key, value in load_2_name_list if value == []]
    for warnId in not_found_list:
        if verbosity > 1:
            log.warning('Following ID has no corresponding entry in the database: %s', warnId)
        else:
            log.debug('Following ID has no corresponding entry in the database: %s', warnId)
    log.info('%s IDs out of %s have not been found', len(not_found_list), len(p_load_list))
    log.info('IDs of missing proteins: %s', not_found_list)
    return not_found_list, load_2_name_list, db_id_list


def erase_custom_fields():
    """
    Resets the custom and main_connex fields of all the Nodes on which we have iterated here.

    Required to be run after node set or node connectivity were modified.
    Unlike the method in the Matrix_retrieval cluster, this method is very time-consuming,
    since it iterates on all the elements of all the classes susceptible to have the custom field.
    """
    def reset_routine(_node):
        _node.custom = ''
        _node.main_connex = False
        _node.save()

    # TODO: redundant implementation - remove one of them
    # TODO: could be accelerated by batching - do it.

    node_gen = DatabaseGraph.find({"custom": "Main_Connex"})
    if len(node_gen) > 0:
        for node in node_gen:
            DatabaseGraph.set_attributes(node.id, {'custom': '', 'main_connex': False})

    node_gen = DatabaseGraph.find({"main_connex": True})
    if len(node_gen) > 0:
        for node in node_gen:
            DatabaseGraph.set_attributes(node.id, {'custom': '', 'main_connex': False})


# TODO: used once elsewhere - consider folding it there
def node_extend_once(edge_type_filter, main_connex_only, core_node):
    """

    :param edge_type_filter:
    :param main_connex_only:
    :param core_node:
    :return:
    """
    node_neighbors = []
    node_neighbor_no = 0
    for edge_type in edge_type_filter:

        for node in DatabaseGraph.get_linked(get_db_id(core_node), link_type=edge_type):
            node_is_connex = node.properties['main_connex']
            if (main_connex_only and node_is_connex) or not main_connex_only:
                node_neo4j_id = get_db_id(node)
                if node_neo4j_id not in forbidden_neo4j_ids:
                    node_neighbors.append(node_neo4j_id)
                    node_neighbor_no += 1

    return node_neighbors, node_neighbor_no


# TODO: add a directionality argument
def expand_from_seed(seed_node_id, edge_filter, main_connex_only):
    """
    Recovers all the nodes accessible in one jump from a seed_node with a given database ID by
    jumping only via the relations of types specified in the edge_filter

    :param seed_node_id: the database ID of the initial node from which we are observing
     accessibility
    :param edge_filter: the list of relation types for which the jumps are authorised
    :param main_connex_only: of true, will expand from seed onto the elements of the main connex
    only
    :return: List of found nodes database ID, number of found nodes
    """
    node_neighbors = []
    for edge_type in edge_filter:
        seed_node_is_connex = DatabaseGraph.get(seed_node_id).properties['main_connex']
        for linked_node in DatabaseGraph.get_linked(seed_node_id, 'both', edge_type):
            if linked_node.id not in forbidden_neo4j_ids and (seed_node_is_connex or not main_connex_only):
                node_neighbors.append(linked_node.id)

    return node_neighbors, len(node_neighbors)

def annotation_ids_from_csv(source_csv):
    """
    Recovers the set of annotation ids from a csv document

    :param source_csv:
    :return:
    """
    analysis_bulbs_ids = []
    with open(source_csv) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            analysis_bulbs_ids = analysis_bulbs_ids + row
    return analysis_bulbs_ids


def cast_analysis_set_to_bulbs_ids(analysis_set_csv_location=analysis_protein_ids_csv):
    """
    Unwraps the bioflow file specified in the analysis_protein_ids_csv, translates its database
    internal ids for further use

    :param analysis_set_csv_location:
    :return:
    """
    analysis_bulbs_ids = annotation_ids_from_csv(analysis_set_csv_location)
    analysis_bulbs_ids = [ret for ret in analysis_bulbs_ids]
    source = look_up_annotation_set(analysis_bulbs_ids)
    PrettyPrinter(indent=4, stream=open(Dumps.analysis_set_display_names, 'w')).pprint(source[1])
    writer(open(Dumps.analysis_set_bulbs_ids, 'w'), delimiter='\n').writerow(source[2])


def cast_background_set_to_bulbs_id(background_set_csv_location=background_protein_ids_csv,
                                    analysis_set_csv_location=analysis_protein_ids_csv):
    """
    Unwraps the bioflow file specified in the background_bulbs_ids, translates it to the database
    internal ids for further use

    :param analysis_set_csv_location:
    :param background_set_csv_location:
    :return:
    """
    background_bulbs_ids = []
    if background_set_csv_location:
        background_bulbs_ids += annotation_ids_from_csv(analysis_set_csv_location)
        background_bulbs_ids += annotation_ids_from_csv(background_set_csv_location)
        background_bulbs_ids = list(set(ret for ret in background_bulbs_ids))

    source = look_up_annotation_set(background_bulbs_ids)
    writer(open(Dumps.background_set_bulbs_ids, 'w'), delimiter='\n').writerow(source[2])


def neo4j_memoize_type(node_type, dict_to_load_into=None):

    if dict_to_load_into is None:
        dict_to_load_into = {}
    log.info('starting %s memoization load', node_type)
    all_nodes_of_type = DatabaseGraph.get_all(node_type)
    for node in all_nodes_of_type:
        dict_to_load_into[node.id] = node
    log.info('%s Loaded, contains %s elements', node_type, len(dict_to_load_into))

    return dict_to_load_into


def recompute_forbidden_ids(forbidden_entities_list):
    """
    Recomputes the list of nodes that contain overloaded terms.

    Without eliminating those terms, they would bring too close together the reactions
    that are normally not, just because of  participation of ultra-abundant elements,
    such as H2O, H+ or ATP

    :param forbidden_entities_list: Dictionary mapping the names of entities to
    their corresponding bulbs classes
    """
    forbidden_ids_list = set()
    for name in forbidden_entities_list:
        for forbidden_legacy_id in Leg_ID_Filter:
            bulbs_class = neo4j_names_dict[name]
            generator = DatabaseGraph.find({"displayName": forbidden_legacy_id}, bulbs_class)

            associated_node_ids = node_generator_2_db_ids(generator)
            forbidden_ids_list.update(associated_node_ids)

    log.info('recomputed %s forbidden IDs. \n Dumping them to %s',
             len(forbidden_ids_list), Dumps.Forbidden_IDs)
    # log.info(forbidden_ids_list)
    pickle.dump(forbidden_ids_list, file(Dumps.Forbidden_IDs, 'w'))


def clear_all(instruction_list):
    """
    empties the whole BioPax-bound node set.

    :param instruction_list:
    """

    for name in instruction_list:
        neo4j_class = neo4j_names_dict[name]
        log.info('processing class: %s', neo4j_class)
        DatabaseGraph.delete_all(neo4j_class)
        log.info('class %s finished processing', neo4j_class)


def run_diagnostics(instructions_list):
    """
    Checks the number of nodes of each type.

    :param instructions_list:
    """
    super_counter = 0
    str_list = ['Database Diagnostics:']
    for name in instructions_list:
        bulbs_class = neo4j_names_dict[name]
        counter = DatabaseGraph.count(bulbs_class)
        str_list.append('\t %s : %s' % (name, counter))
        super_counter += counter
    str_list.append('Total : %s' % super_counter)
    log.info('\n'.join(str_list))


def convert_to_internal_ids(base):
    """"
    converts names of proteins to database_ids, preferably matching to UNIPROTS

    :param base:
    :return:
    """

    warn_list, results_tuple_list, results_list = look_up_annotation_set(set(base))
    return_dict = {}

    breakpoints = 300
    size = len(results_tuple_list)

    for i, (key, match_list) in enumerate(results_tuple_list):
        if i % breakpoints == 0:
            log.info("\t %.2f %%" % (float(i)/float(size)*100))
        if key not in warn_list:
            for match in match_list:
                if match[0] == 'UNIPROT':
                    return_dict[key] = match[2]
                else:
                    return_dict[key] = match_list[0][2]

    log.debug('ID cast converter length: %s', len(return_dict))

    return return_dict


def cross_link_identifiers():
    log.info('Cross-linking the identifiers: Uniprot Acnums')
    DatabaseGraph.cross_link_on_annotations('UNIPROT_Accnum')
    log.info('Cross-linking the identifiers: Uniprot Names')
    DatabaseGraph.cross_link_on_annotations('UNIPROT_Name')
    log.info('Cross-linking the identifiers: Uniprot Gene Names')
    DatabaseGraph.cross_link_on_annotations('UNIPROT_GeneName')


# Yes, I know what goes below here is ugly and shouldn't be in the
# production part of the code

on_rtd = os.environ.get('READTHEDOCS') == 'True'
on_unittest = os.environ.get('UNITTESTING') == 'True'

if on_unittest:
    # Yes, this is dangerous as hell. Can't see a better way of doing it
    # though for now.
    import sys
    import unittests.Mocks.DB_IO_Mocks as SelfMock
    sys.modules[__name__] = SelfMock


if __name__ == "__main__":
    # erase_custom_fields()
    # recompute_forbidden_ids(forbidden_verification_list)
    # print cast_analysis_set_to_bulbs_ids()
    # print look_up_annotation_set('CTR86')
    # print look_up_annotation_set('ENSG00000131981', 'UNIPROT_Ensembl')
    # cast_analysis_set_to_bulbs_ids()
    # cast_background_set_to_bulbs_id()
    _, resdict, reslist = look_up_annotation_set(['RNF14'])
    pprint(resdict)
    print reslist
    _, resdict, reslist = look_up_annotation_set(['RNF14'])
    pprint(resdict)
    print reslist
    # run_diagnostics(full_list)
    # memoize_bulbs_type(neo4j_names_dict['UNIPROT'][0])
    # cast_analysis_set_to_bulbs_ids()
    # cast_background_set_to_bulbs_id(background_set_csv_location=None)
    # Akshay_p53_go_set = ["0030330", "0000019", "0000002", "0006977"]
    # go_bulbs_ids = lookup("GO Term", Akshay_p53_go_set)
    # association_dict = pull_associations(go_bulbs_ids, ["is_go_annotation"])
    # proteins = [item for key, val in association_dict.iteritems() for item in val]
    # writer(open(Dumps.analysis_set_bulbs_ids, 'w'), delimiter='\n').writerow(proteins)
    # print proteins

    # bulbs_id = lookup('UNIPROT', 'Q16206')
