"""
Responsible for the injection of the transcription factor databases into the main space
"""

from bioflow.bio_db_parsers.tfParsers import parse_marbach, parse_cellnet_grn, parse_TRRUST
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import marbach_mode, marbach_path, cellnet_path, trrust_path
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
    log.debug('TF ID cast converter length: %s', len(return_dict))
    return return_dict


def insert_into_the database()