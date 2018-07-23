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
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph, on_alternative_graph
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
    if on_alternative_graph:
        return neo4j_node.id
    else:
        return neo4j_node._id


# TODO: REFACTOR: offload annotations to ElasticSearch engine
def get_attached_annotations(neo4j_node_id):
    """
    Recovers ids of annotation nodes attached the the node with a given bulbs id

    :param neo4j_node_id:
    :return:
    """
    list_of_annotations = []
    if on_alternative_graph:
        annotation_node_generator = DatabaseGraph.get_linked(neo4j_node_id, link_type="is_annotated")
    else:
        node = DatabaseGraph.vertices.get(neo4j_node_id)
        annotation_node_generator = node.outV("is_annotated")
    if not annotation_node_generator:
        log.debug("node %s has not annotations attached to it", neo4j_node_id)
        return []
    else:
        for rel_node in annotation_node_generator:
            list_of_annotations.append(get_db_id(rel_node))
        return list_of_annotations

# used in bulbs exclusively
def annotations_2_node_chars(node_generator):
    """
    Gets nodes that are referenced by the annot_nodes in provided generator.

    Returns the characterisation of nodes referenced from the node generator

    :param node_generator: iterator over annot_nodes
    :return: the list of object nodes accessible from this set of annot_nodes. [type, name,
    bulbs_id, legacy_id]
    """

    if not node_generator:  # Node generator is empty: return empty list
        return []

    referenced_nodes = []
    for node in node_generator:
        annotates_backlink_generator = node.inV("is_annotated")

        if not annotates_backlink_generator:
            log.warning("%s is floating alone in the wild. He feels lonely. %s",
                        str(node), 'The database is most likely broken.')

        else:
            for object_node in annotates_backlink_generator:
                node_bulbs_id = get_db_id(object_node)
                node_legacy_id = object_node.ID
                node_type = object_node.element_type
                node_display_name = object_node.displayName
                referenced_nodes.append((node_type, node_display_name,
                                         node_bulbs_id, node_legacy_id))

    return referenced_nodes


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

    def double_index_search(pay_load, payload_type):
        """
        Supporting function. Performs a search in the annotation nodes over both a payload and a
        payload type

        :param pay_load: payload
        :param payload_type: payload type
        :return: list of found annotation nodes satisfying both payload content and payload type
         conditions
        """
        inner_node_generator = DatabaseGraph.AnnotNode.index.lookup(payload=pay_load)
        results = []
        if inner_node_generator:  # if there is anything in the inner node generator
            for _node in inner_node_generator:
                if payload_type in _node.ptype:  # late filtering by the node ptype
                    results.append(_node)
        return results

    payload = p_load.upper()  # make payload case-agnostic

    if on_alternative_graph:
        nodes = DatabaseGraph.get_from_annotation_tag(payload, p_type)
        retlist = []
        for node in nodes:
            node_bulbs_id = get_db_id(node)
            node_legacy_id = node.properties['legacyId']
            node_type = list(node.labels)[0]
            node_display_name = node.properties['displayName']
            retlist.append((node_type, node_display_name, node_bulbs_id, node_legacy_id))
        return retlist

    else:
        if p_type in annotation_nodes_ptypes:
            node_generator = double_index_search(payload, p_type)
            return annotations_2_node_chars(node_generator)

        if p_type != '':
            log.exception("%s is an unsupported payload type. ", p_type)
            log.info("Please refer authorised annotation_nodes_ptypes in cofigs/internal_configs.py.")
            log.info("ignoring p_type from now on")

        node_generator = DatabaseGraph.AnnotNode.index.lookup(payload=payload)
        return annotations_2_node_chars(node_generator)


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

    if on_alternative_graph:  # In the end, it was not worth it...
        load_2_name_list = [(p_load, transform_annotation_node(p_nodes)) for (p_load, p_nodes) in
                            zip(p_load_list, DatabaseGraph.batch_retrieve_from_annotation_tags(p_load_list, p_type))]
    else:
        load_2_name_list = [(p_load, look_up_annotation_node(p_load, p_type))
                            for p_load in p_load_list]

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

    if on_alternative_graph:
        node_gen = DatabaseGraph.find({"custom": "Main_Connex"})
        if len(node_gen) > 0:
            for node in node_gen:
                DatabaseGraph.set_attributes(node.id, {'custom': '', 'main_connex': False})

        node_gen = DatabaseGraph.find({"main_connex": True})
        if len(node_gen) > 0:
            for node in node_gen:
                DatabaseGraph.set_attributes(node.id, {'custom': '', 'main_connex': False})

    else:
        node_gen = DatabaseGraph.vertices.index.lookup(costum='Main_Connex')
        if node_gen:
            for node in node_gen:
                reset_routine(node)

        node_gen = DatabaseGraph.vertices.index.lookup(main_connex=True)
        if node_gen:
            for node in node_gen:
                reset_routine(node)


# TODO: delete: this should be dead code in the final refactoring version
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

        if on_alternative_graph:
            for node in DatabaseGraph.get_linked(get_db_id(core_node), link_type=edge_type):
                node_is_connex = node.properties['main_connex']
                if (main_connex_only and node_is_connex) or not main_connex_only:
                    node_neo4j_id = get_db_id(node)
                    if node_neo4j_id not in forbidden_neo4j_ids:
                        node_neighbors.append(node_neo4j_id)
                        node_neighbor_no += 1

        else:
            if core_node.bothV(edge_type) is not None:  # get the next edge type

                # TODO: REFACTOR: connex logic is quite unclear
                for node in core_node.bothV(edge_type):  # get all the nodes in core_node
                    connex = True   # by default all nodes are connex

                    if main_connex_only:  # when we are looking for main connex it is not in the main connex
                        connex = False

                        if node.main_connex:  # unless the node from which we are connecting is connex
                            connex = True

                    node_bulbs_id = get_db_id(node)  # get bulbs ID of a new node

                    if node_bulbs_id not in forbidden_neo4j_ids and connex:
                        node_neighbors.append(node_bulbs_id)
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
    if on_alternative_graph:
        node_neighbors = []
        for edge_type in edge_filter:
            seed_node_is_connex = DatabaseGraph.get(seed_node_id).properties['main_connex']
            for linked_node in DatabaseGraph.get_linked(seed_node_id, 'both', edge_type):
                if linked_node.id not in forbidden_neo4j_ids and (seed_node_is_connex or not main_connex_only):
                    node_neighbors.append(linked_node.id)

        return node_neighbors, len(node_neighbors)

    else:
        seed_node = DatabaseGraph.vertices.get(seed_node_id)
        return node_extend_once(edge_filter, main_connex_only, seed_node)


#TODO : dead code, delete it
def _recover_annotation(node_id_set, annotation_type):
    """
    seems to retrieve annotations attached to a list of nodes in a set.

    :param node_id_set:
    :param annotation_type:
    :return:
    """

    ret_dict = defaultdict(list)
    ret_list = []

    if on_alternative_graph:
        for node_id in node_id_set:
            annotation_nodes = DatabaseGraph.get_linked(node_id, 'to', 'is_annotated')
            for annotation_node in annotation_nodes:
                if annotation_node.properties["type"] == annotation_type:
                    ret_dict[node_id].append(annotation_node.properties["tag"])
                    ret_list.append((annotation_node.properties["tag"], annotation_node.id))

    else:
        for node_id in node_id_set:

            node = DatabaseGraph.vertices.get(node_id)

            if node:
                annotation_generator = node.bothV('is_annotated')

                if annotation_generator:
                    for annot in annotation_generator:

                        if annot.ptype == annotation_type:
                            if 'G0' in annot.payload:  # What? Why?
                                ret_dict[node_id].append(annot.payload)
                                ret_list.append([str(annot.payload), str(node.ID)])

    return dict(ret_dict), ret_list


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


# TODO: bulbflow-specific workaround; delete it
def _bulb_specific_stable_get_all(bulbs_class):
    """
    A Bulbs_NodeProxy.get_all() wrapper that guards against the failure of the neo4j index.

    :param bulbs_class: Node Proxy object
    """
    unsafe_generator = bulbs_class.get_all()
    type_var = bulbs_class.client.config.type_var
    element_type = bulbs_class.element_class.get_element_type(bulbs_class.client.config)
    if unsafe_generator:
        for bulbs_node in unsafe_generator:
            if getattr(bulbs_node, type_var) == element_type:
                yield bulbs_node


# TODO: bulbflow-specific implementation; delete it
def _bulb_specific_memoize_bulbs_type(bulbs_class, dict_to_load_into=None):
    """
    Loads a bulbs type ID_2_object into a supplied dict. If no dict supplied, returns the dict as
    a result

    :param bulbs_class: bulbs class whose contents we want to memoize
    :param dict_to_load_into: if provided, will be loaded into and then returned
    :return:
    """
    class_name = bulbs_class.element_class.get_element_type(bulbs_class.client.config)
    if dict_to_load_into is None:
        dict_to_load_into = {}
    log.info('starting %s memoization load', class_name)

    for bulbs_node in _bulb_specific_stable_get_all(bulbs_class):
        dict_to_load_into[bulbs_node.ID] = bulbs_class.get(get_db_id(bulbs_node))
    log.info('%s Loaded, contains %s elements', class_name, len(dict_to_load_into))

    return dict_to_load_into


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
        bulbs_class, _ = neo4j_names_dict[name]
        for forbidden_legacy_id in Leg_ID_Filter:
            if on_alternative_graph:
                generator = DatabaseGraph.find({"displayName": forbidden_legacy_id}, bulbs_class)
            else:
                generator = bulbs_class.index.lookup(displayName=forbidden_legacy_id)

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
    def bulbs_clear_all():
        for name in instruction_list:
            bulbs_class, bulbs_alias = neo4j_names_dict[name]
            log.info('processing class: %s, alias %s', bulbs_class, bulbs_alias)
            if _bulb_specific_stable_get_all(bulbs_class):
                id_list = [get_db_id(bulbs_class_instance)
                           for bulbs_class_instance in _bulb_specific_stable_get_all(bulbs_class)]
                id_list_length = len(id_list) / 100.
                for counter, ID in enumerate(id_list):
                    del_set = get_attached_annotations(ID)
                    for annotation_node_id in del_set:
                        DatabaseGraph.AnnotNode.delete(annotation_node_id)
                    bulbs_class.delete(ID)
                    if counter % 100 == 0:
                        log.info('deleting class %s %.2f %%:', name, counter / id_list_length)
                log.info('deleting class %s %.2f %%:', name, 100)

    def cypher_clear_all():
        for name in instruction_list:
            neo4j_class = neo4j_names_dict[name]
            log.info('processing class: %s', neo4j_class)
            DatabaseGraph.delete_all(neo4j_class)
            log.info('class %s finished processing', neo4j_class)

    if on_alternative_graph:
        cypher_clear_all()
    else:
        bulbs_clear_all()


def run_diagnostics(instructions_list):
    """
    Checks the number of nodes of each type.

    :param instructions_list:
    """
    def bulbs_diagnostics():
        super_counter = 0
        str_list = ['Database Diagnostics:']
        for name in instructions_list:
            bulbs_class, bulbs_alias = neo4j_names_dict[name]
            counter = bulbs_class.index.count(element_type=bulbs_alias)
            str_list.append('\t %s : %s' % (name, counter))
            super_counter += counter
        str_list.append('Total : %s' % super_counter)
        log.info('\n'.join(str_list))

    def cypher_diagnostics():
        super_counter = 0
        str_list = ['Database Diagnostics:']
        for name in instructions_list:
            bulbs_class = neo4j_names_dict[name]
            counter = DatabaseGraph.count(bulbs_class)
            str_list.append('\t %s : %s' % (name, counter))
            super_counter += counter
        str_list.append('Total : %s' % super_counter)
        log.info('\n'.join(str_list))

    if on_alternative_graph:
        cypher_diagnostics()
    else:
        bulbs_diagnostics()

#TODO: dead code, just delete it
def lookup(bulbs_node_shortname, list_of_ids):
    """
    Just performs the lookup for a list of legacy IDs in for a given node type

    :param bulbs_node_shortname:
    :param list_of_ids:
    :return:
    """
    def bulbs_lookup():
        bulbs_node_id_list = []
        for legacy_id in list_of_ids:
            bulbs_node_generator = neo4j_names_dict[bulbs_node_shortname][0].index.lookup(ID=legacy_id)
            if bulbs_node_generator:
                for node in bulbs_node_generator:
                    bulbs_node_id_list.append(get_db_id(node))
                    print node.ID, get_db_id(node), node.displayName
        return bulbs_node_id_list

    def cypher_lookup():
        neo4j_nodes_id_list = []
        for legacy_id in list_of_ids:
            for node in DatabaseGraph.find({'legacyID': legacy_id}):
                neo4j_nodes_id_list.append(get_db_id(node))
                print node.legacyID, get_db_id(node), node.parameters['displayName']
        return neo4j_nodes_id_list

    if on_alternative_graph:
        return cypher_lookup()
    else:
        return bulbs_lookup()

# TODO: dead code, just detelete it
def pull_associations(bulbs_nodes_id_list, allowed_edge_types, direction="BOTH"):
    """
    Pulls associated nodes to a given set of bulbs IDs according to a given relationship

    :param bulbs_nodes_id_list:
    :param allowed_edge_types:
    :param direction:
    :return:
    """
    def bulbs_pull_associations():
        association_dict = {}
        for bulbs_node_id in bulbs_nodes_id_list:
            root_node = DatabaseGraph.vertices.get(bulbs_node_id)
            association_dict[bulbs_node_id] = node_extend_once(allowed_edge_types,
                                                               False, root_node)[0]
        return association_dict

    def cypher_pull_associations():
        association_dict = {}
        for bulbs_node_id in bulbs_nodes_id_list:
            association_dict[bulbs_node_id] = node_extend_once(allowed_edge_types, False, bulbs_node_id)

    if on_alternative_graph:
        return cypher_pull_associations()
    else:
        return bulbs_pull_associations()


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
