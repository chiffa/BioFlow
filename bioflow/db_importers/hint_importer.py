"""
Set of tools to work with HiNT database
"""
from bioflow.bio_db_parsers.proteinRelParsers import parse_hint
from bioflow.main_configs import hint_csv_path
from bioflow.neo4j_db.db_io_routines import _bulb_specific_memoize_bulbs_type
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def get_uniprots_for_hint():
    """
    Recovers UP Gene names maps to UNIPROT nodes containing them.

    :return:
    """
    inital_dict = _bulb_specific_memoize_bulbs_type(DatabaseGraph.UNIPORT)
    for key in inital_dict.keys():
        inital_dict[key.split('_')[0]] = inital_dict.pop(key)
    return inital_dict


def cross_ref_hint(flush=True):
    """
    Pulls Hint relationships and connects reached_uniprots_bulbs_id_list in the database

    :param flush: if True, relationships are pushed to the actual graph database
    :return:
    """
    relations_dict = parse_hint(hint_csv_path)
    uniprot_ref_dict = get_uniprots_for_hint()
    processed_pairs = set()
    actual_cross_links = 0
    for key in uniprot_ref_dict.keys():
        if key in relations_dict.keys():
            processed_pairs.add(key)
            for subkey in relations_dict[key]:
                if subkey in uniprot_ref_dict.keys() and subkey not in processed_pairs:
                    actual_cross_links += 1
                    log.debug('HINT links: %s, %s', key, subkey)
                    if flush:
                        DatabaseGraph.is_interacting.create(
                            uniprot_ref_dict[key], uniprot_ref_dict[subkey], source='HINT')
    log.info('HINT Cross-links: %s, HINT processed pairs: %s',
             actual_cross_links, len(processed_pairs))
