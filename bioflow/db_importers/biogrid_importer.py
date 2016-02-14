"""
Responsible for the injection of the BioGRID parse into the main space.
"""
from bioflow.bio_db_parsers.proteinRelParsers import parse_bio_grid
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import biogrid_path
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph


log = get_logger(__name__)


def convert_to_internal_ids(base):
    """
    converts names of proteins to database ids

    :param base:
    :return:
    """
    warn_list, results_tuple_list, results_list = look_up_annotation_set(set(base))
    return_dict = dict((key, value[0][2])
                       for key, value in results_tuple_list if key not in warn_list)
    log.debug('BioGrid ID cast converter length: %s', len(return_dict))
    return return_dict


def insert_into_the_database(_up_ids_2_inner_ids, _up_ids_2_properties):
    """
    performs the insertion into the database subroutine

    :param _up_ids_2_inner_ids:
    :param _up_ids_2_properties:
    :return:
    """

    final_dicts = dict(
        ((_up_ids_2_inner_ids[key[0]], _up_ids_2_inner_ids[key[1]]), value)
        for key, value in _up_ids_2_properties.iteritems()
        if key[0] in _up_ids_2_inner_ids.keys() and key[1] in _up_ids_2_inner_ids.keys())

    for (node1_id, node2_id), link_parameters in final_dicts.iteritems():
        node1 = DatabaseGraph.UNIPORT.get(node1_id)
        node2 = DatabaseGraph.UNIPORT.get(node2_id)
        if len(link_parameters) > 1:
            DatabaseGraph.is_weakly_interacting.create(node1, node2,
                                                       throughput=link_parameters[0],
                                                       confidence=float(link_parameters[1]))
        else:
            DatabaseGraph.is_weakly_interacting.create(node1, node2,
                                                       throughput=link_parameters[0])


def cross_ref_bio_grid():
    """    performs the total biogrid_path parse. """
    log.info('Starting BioGrid Parsing')
    up_ids_2_properties, up_ids = parse_bio_grid(biogrid_path)
    log.info('BioGrid parsed, starting translation of UP identifiers to internal database ' +
             'identifiers')
    up_ids_2_inner_ids = convert_to_internal_ids(up_ids)
    log.info('UP identifier conversion finished, starting database insertion')
    insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties)
    log.info('Database insertion finished. biogrid_path import finished')
