"""
Responsible for the injection of the BioGRID parse into the main space.
"""
from bioflow.bio_db_parsers.proteinRelParsers import parse_bio_grid
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import biogrid_path
from bioflow.neo4j_db.db_io_routines import convert_to_internal_ids
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from time import time
import datetime


log = get_logger(__name__)


def insert_into_the_database(_up_ids_2_inner_ids, _up_ids_2_properties):
    """
    performs the insertion into the database subroutine

    :param _up_ids_2_inner_ids:
    :param _up_ids_2_properties:
    :return:
    """

    final_dicts = dict(
        ((_up_ids_2_inner_ids[key[0]], _up_ids_2_inner_ids[key[1]]), value)
        for key, value in _up_ids_2_properties.items()
        if key[0] in list(_up_ids_2_inner_ids.keys()) and key[1] in list(_up_ids_2_inner_ids.keys()))

    breakpoints = 300
    total_pairs = len(list(final_dicts.keys()))
    previous_time = time()

    for counter, ((node1_id, node2_id), link_parameters) in enumerate(final_dicts.items()):

        if counter % breakpoints == 0 and counter > 1:
            compops = float(breakpoints)/(time() - previous_time)
            mins_before_termination = int((total_pairs-counter)/compops/60)

            log.info("\t %.2f %% done, finishing in %s min" % (float(counter)/float(total_pairs)*100, mins_before_termination))

            log.debug('inserting link %s out of %s; %.2f complete; inserting speed: %.2f; expected finsihing: %s mins',
                      counter + 1,
                      total_pairs,
                      counter / float(total_pairs) * 100,
                      compops,
                      datetime.datetime.now() + datetime.timedelta(minutes=mins_before_termination))
            previous_time = time()

        if len(link_parameters) > 1:
            DatabaseGraph.link(node1_id, node2_id, 'is_weakly_interacting',
                               {'source': 'BioGRID', 'throughput': link_parameters[0],
                                'confidence': float(link_parameters[1])})
        else:
            DatabaseGraph.link(node1_id, node2_id, 'is_weakly_interacting',
                               {'source': 'BioGRID', 'throughput': link_parameters[0]})



def cross_ref_bio_grid():
    """    performs the total BioGRID parse. """
    log.info('Starting BioGRID Parsing')
    up_ids_2_properties, up_ids = parse_bio_grid(biogrid_path)
    log.info('BioGrid parsed, starting translation of UP identifiers to internal database ' +
             'identifiers')
    up_ids_2_inner_ids = convert_to_internal_ids(up_ids)
    log.info('UP identifier conversion finished, starting database insertion')
    insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties)
    log.info('Database insertion finished. BioGRID import finished')
