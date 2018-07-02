"""
Imports the complexes and insert them into the neo4j database
"""
from bioflow.bio_db_parsers.ComplexPortalParser import parse_complex_portal
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import complexes_path
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from time import time
import numpy as np
import datetime


log = get_logger(__name__)


def convert_to_internal_ids(base):
    """"
    converts names of proteins to database_ids

    :param base:
    :return:
    """
    warn_list, results_tuple_list, results_list = look_up_annotation_set(set(base))
    return_dict = dict((key, value[0][2])
                       for key, value in results_tuple_list if key not in warn_list)
    log.debug('TF ID cast converter length: %s', len(return_dict))
    return return_dict


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
    for counter, new_node in enumerate(up_ids_2_properties.iteritems()):

        complex_node = DatabaseGraph.COMPLEX.create(ID=new_node['ID'],
                                                    displayName=new_node['DisplayName'],
                                                    main_connex=False)

        for node2_up in new_node['components']:

            node2 = DatabaseGraph.UNIPORT.get(up_ids_2_inner_ids[node2_up])

            DatabaseGraph.is_interacting.create(complex_node, node2,
                                                source=origin,
                                                weight=1.0)


    final_dicts = dict(
        ((up_ids_2_inner_ids[key[0]], up_ids_2_inner_ids[key[1]]), value)
        for key, value in up_ids_2_properties.iteritems()
        if key[0] in up_ids_2_inner_ids.keys() and key[1] in up_ids_2_inner_ids.keys())

    breakpoints = 300
    total_pairs = len(final_dicts.keys())
    previous_time = time()

    for counter, ((node1_id, node2_id), link_parameter) in enumerate(final_dicts.iteritems()):
        # print 'tick'

        if counter % breakpoints == 0 and counter > 1:
            compops = float(breakpoints)/(time()-previous_time)
            secs_before_termination = int((total_pairs-counter)/compops)

            log.info('inserting link %s out of %s; %.2f complete; inserting speed: %.2f; expected finsihing: %s',
                      counter + 1,
                      total_pairs,
                      counter / float(total_pairs) * 100,
                      compops,
                      datetime.datetime.now() + datetime.timedelta(seconds=secs_before_termination))
            previous_time = time()

        node1 = DatabaseGraph.UNIPORT.get(node1_id)
        node2 = DatabaseGraph.UNIPORT.get(node2_id)

        DatabaseGraph.is_interacting.create(node1, node2,
                                            source=origin,
                                            weight=float(np.sum(np.array(link_parameter).astype(np.int))))


def insert_complexes():
    """
    Performs the full kinase-substrate parsing and insertion.

    :return:
    """
    log.info('Starting Complex Portal parsing')

    up_ids_2_properties, up_ids = parse_complex_portal(complexes_path)

    log.info('Complex Portal parsed, starting translation of UP identifiers to internal database identifiers')
    up_ids_2_inner_ids = convert_to_internal_ids(up_ids)

    log.info('UP identifier conversion finished, starting database insertion for %s links' % len(up_ids_2_properties.keys()))
    insert_into_the_database(up_ids_2_inner_ids, up_ids_2_properties, 'ComplexPortal')

    log.info('Database insertion finished. Complex Portal import finished')