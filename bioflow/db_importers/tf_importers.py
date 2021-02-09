"""
Responsible for the injection of the transcription factor databases into the main space
"""

from bioflow.bio_db_parsers.tfParsers import parse_TRRUST
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import trrust_path, organism
from bioflow.neo4j_db.db_io_routines import convert_to_internal_ids
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from time import time
import datetime


log = get_logger(__name__)


def insert_into_the_database(up_ids_2_inner_ids,
                             up_ids_2_properties,
                             strong_interaction_threshold,
                             origin):
    """
    Performs the insertion in the database sub-routine

    :param up_ids_2_inner_ids:
    :param up_ids_2_properties:
    :param strong_interaction_threshold:
    :param origin:
    :return:
    """

    final_dicts = dict(
        ((up_ids_2_inner_ids[key[0]], up_ids_2_inner_ids[key[1]]), value)
        for key, value in up_ids_2_properties.items()
        if key[0] in list(up_ids_2_inner_ids.keys()) and key[1] in list(up_ids_2_inner_ids.keys()))

    breakpoints = 300
    total_pairs = len(list(final_dicts.keys()))
    previous_time = time()

    for counter, ((node1_id, node2_id), link_parameter) in enumerate(final_dicts.items()):
        # print 'tick'

        if counter % breakpoints == 0 and counter > 1:
            compops = float(breakpoints) / (time() - previous_time)
            secs_before_termination = int((total_pairs - counter) / compops)

            log.info('inserting link %s out of %s; %.2f complete; inserting speed: %.2f; expected finsihing: %s',
                      counter + 1,
                      total_pairs,
                      counter / float(total_pairs) * 100,
                      compops,
                      datetime.datetime.now() + datetime.timedelta(seconds=secs_before_termination))
            previous_time = time()

        if link_parameter >= strong_interaction_threshold:
            DatabaseGraph.link(node1_id, node2_id, 'is_interacting',
                               {'source': origin,
                                'weight': float(link_parameter)})
        else:
            DatabaseGraph.link(node1_id, node2_id, 'is_weakly_interacting',
                               {'source': origin,
                                'weight': float(link_parameter)})


def cross_ref_tf_factors(confs='tcm'):
    """
    Performs the full transcription factors parsing and insertion routine.

    :return:
    """
    # TODO:
    if organism != 'Human':
        raise Exception('TF data unavailable for organisms other than human. Disable TF import in import_main.py')

    if 't' in confs:
        log.info('Starting TRRUST Parsing')
        up_ids_2_properties, up_ids = parse_TRRUST(trrust_path)

        log.info('TRRUST parsed, starting translation of UP identifiers to internal database identifiers')
        up_ids_2_inner_ids = convert_to_internal_ids(up_ids)

        log.info('UP identifier conversion finished, starting database insertion for %s links' % len(list(up_ids_2_properties.keys())))
        insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties, 1, 'TF_TRRUST')

        log.info('Database insertion finished. TRRUST import finished')

    # if 'c' in confs:
    #     log.info('Starting CellNet Parsing')
    #     up_ids_2_properties, up_ids = parse_cellnet_grn(cellnet_path)
    #
    #     log.info('CellNet parsed, starting translation of UP identifiers to internal database identifiers')
    #     up_ids_2_inner_ids = convert_to_internal_ids(up_ids)
    #
    #     log.info('UP identifier conversion finished, starting database insertion for %s links' % len(list(up_ids_2_properties.keys())))
    #     insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties, 10, 'CellNet')
    #
    #     log.info('Database insertion finished. CellNet import finished')

    # if 'm' in confs:
    #     log.info('Starting Marbach Parsing')
    #     up_ids_2_properties, up_ids = parse_marbach(marbach_path, marbach_mode)
    #
    #     log.info('Marbach parsed, starting translation of UP identifiers to internal database identifiers')
    #     up_ids_2_inner_ids = convert_to_internal_ids(up_ids)
    #
    #     log.info('UP identifier conversion finished, starting database insertion for %s links' % len(list(up_ids_2_properties.keys())))
    #     insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties, 10, 'Marbach2016')
    #
    #     log.info('Database insertion finished. Marbach import finished')
