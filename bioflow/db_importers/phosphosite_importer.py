"""
Responsible for the injection of the phosphosite phosphorilation regulation into the main space
"""
from bioflow.bio_db_parsers.PhosphositeParser import parse_phosphosite
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import phosphosite_path, phosphosite_organism
from bioflow.neo4j_db.db_io_routines import convert_to_internal_ids
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from time import time
import numpy as np
import datetime


log = get_logger(__name__)


def insert_into_the_database(up_ids_2_inner_ids,
                             up_ids_2_properties,
                             origin):
    """
    Performs the insertion in the database sub-routine

    :param up_ids_2_inner_ids:
    :param up_ids_2_properties:
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

        if counter % breakpoints == 0 and counter > 1:
            # TODO: [progress bar]

            compops = float(breakpoints) / (time() - previous_time)
            secs_before_termination = int((total_pairs - counter) / compops)

            log.info('inserting link %s out of %s; %.2f complete; inserting speed: %.2f; expected finsihing: %s',
                      counter + 1,
                      total_pairs,
                      counter / float(total_pairs) * 100,
                      compops,
                      datetime.datetime.now() + datetime.timedelta(seconds=secs_before_termination))
            previous_time = time()

        DatabaseGraph.link(node1_id, node2_id,
                           'is_interacting',
                           {'source': origin,
                            'weight': float(np.sum(np.array(link_parameter).astype(int))),
                            'parse_type': 'physical_entity_molecular_interaction'})


def cross_ref_kinases_factors():
    """
    Performs the full kinase-substrate parsing and insertion.

    :return:
    """
    log.info('Starting PhosphoSite Parsing')

    up_ids_2_properties, up_ids = parse_phosphosite(phosphosite_path, phosphosite_organism)

    log.info('PhosphoSite parsed, starting translation of UP identifiers to internal database identifiers')
    up_ids_2_inner_ids = convert_to_internal_ids(up_ids)

    log.info('UP identifier conversion finished, starting database insertion for %s links' % len(list(up_ids_2_properties.keys())))
    insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties, 'PhosphoSite')

    log.info('Database insertion finished. PhosphoSite import finished')
