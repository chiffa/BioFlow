# TODO: functions from here are implemented elsewhere: module to be deleted

from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph, on_alternative_graph
from bioflow.utils.log_behavior import get_logger
from bioflow.neo4j_db.db_io_routines import _bulb_specific_stable_get_all

log = get_logger(__name__)


def print_node_look_up_by_id(domain, search_legacy_id):
    """
    Looks up a node by legacy ID and prints it's access url. The lookup is a case-sensitive strict
    match.

    :param domain: Node type of the form of DatabaseGraph.Object
    :param search_legacy_id: requested legacy ID
    :return: list of the items found
    """
    def bulbs_way():
        accumulator = []
        lookup_generator = domain.index.lookup(ID=search_legacy_id)
        if not lookup_generator:
            print "nothing found"
        if lookup_generator:
            for item in lookup_generator:
                print item
                accumulator.append(item)
        return accumulator

    def native_drive_way():
        accumulator = []
        for item in DatabaseGraph.find({'legacyId':search_legacy_id}):
            print item
            accumulator.append(item)
        return accumulator

    if on_alternative_graph:
        return native_drive_way()
    else:
        return bulbs_way()


def count_nodes(domain):
    """
    Stupid counter that gets the number of items in a given domain

    :param domain: Domain of the DatabaseGraph.Object whose objects we are going to count
    :return: number of objects in the domain
    """
    if on_alternative_graph:
        return DatabaseGraph.count(domain)
    else:
        return sum(1 for _ in _bulb_specific_stable_get_all(domain))


def print_look_up_by_id_for_a_set(domain, id_set):
    """
    Looks up nodes by legacy IDs from an ID_set list and prints their access url. The lookup is a
    case-sensitive strict match.

    :param domain: Node type of the form of DatabaseGraph.Object
    :param id_set: list of requested legacy ID
    """
    for ID in id_set:
        print "scanning for:", ID
        print_node_look_up_by_id(domain, ID)
        print "=================="


def get_display_names(up_nodes, are_neo4j_ids):
    """
    Maps Uniprot nodes legacy or bulbs IDs to their display names

    :param up_nodes:
    :param are_neo4j_ids:
    :return:
    """
    up_legacy_ids_2_display_names = {}

    if are_neo4j_ids:
        for node_neo4j_id in up_nodes:
            if on_alternative_graph:
                node = DatabaseGraph.get(node_neo4j_id)
                up_legacy_ids_2_display_names[node.properties['legacyId']] = [node.properties['displayName']]

            else:
                node = DatabaseGraph.UNIPORT.get(node_neo4j_id)
                up_legacy_ids_2_display_names[node.ID] = [node.displayName]
        return up_legacy_ids_2_display_names

    for node_legacy_id in up_nodes:
        if on_alternative_graph:
            up_legacy_ids_2_display_names[node_legacy_id] = []
            for node in DatabaseGraph.find({'legacyId': node_legacy_id}):
                up_legacy_ids_2_display_names[node_legacy_id].append(node.properties['displayName'])

        else:
            generator = DatabaseGraph.UNIPORT.index.lookup(ID=node_legacy_id)
            if not generator:
                continue
            display_names_list = []
            for node in generator:
                display_names_list.append(node.displayName)
            if len(display_names_list) != 1:
                log.exception('Something wrong with UP %s; too many display names: %s',
                              node_legacy_id, display_names_list)
            else:
                up_legacy_ids_2_display_names[node_legacy_id] = display_names_list

    return up_legacy_ids_2_display_names


if __name__ == "__main__":
    print count_nodes(DatabaseGraph.UNIPORT)
    print_node_look_up_by_id(DatabaseGraph.UNIPORT, "SIR2_YEAST")
