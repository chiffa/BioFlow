"""
Responsible for the injection of the BioGRID parse into the main space.
"""
from BioFlow.bio_db_parsers.proteinRelParsers import parse_bio_grid
from BioFlow.utils.log_behavior import logger
from BioFlow.main_configs import BioGRID
from BioFlow.neo4j_db.db_io_routines import look_up_annotation_set
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph


def convert_to_internal_ids(base):
    """
    converts names of proteins to database ids

    :param base:
    :return:
    """
    warn_list, results_dict, results_list = look_up_annotation_set(set(base))
    return_dict = dict((key, value[0][2])
                       for key, value in results_dict.iteritems() if key not in warn_list)
    logger.debug('BioGrid ID cast converter length: %s' % len(return_dict))
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
        # TODO: This can be factored out to the database io as cross-linking Uniprots
        node1 = DatabaseGraph.UNIPORT.get(node1_id)
        node2 = DatabaseGraph.UNIPORT.get(node2_id)
        if len(link_parameters) > 1:
            DatabaseGraph.is_weakly_interacting.create(node1, node2,
                                                       throughput=link_parameters[0],
                                                       confidence=float(link_parameters[1]))
        else:
            DatabaseGraph.is_weakly_interacting.create(node1, node2,
                                                       throughput=link_parameters[0])


def import_bio_grid():
    """    performs the total BioGRID parse. """
    logger.info('Starting BioGrid Parsing')
    up_ids_2_properties, up_ids = parse_bio_grid(BioGRID)
    logger.info('BioGrid parsed, starting translation of UP identifiers to internal database ' +
                'identifiers')
    up_ids_2_inner_ids = convert_to_internal_ids(up_ids)
    logger.info('UP identifier conversion finished, starting database insertion')
    insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties)
    logger.info('Database insertion finished. BioGRID import finished')
