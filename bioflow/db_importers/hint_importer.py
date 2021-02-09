"""
Set of tools to work with HiNT database
"""
from bioflow.bio_db_parsers.proteinRelParsers import parse_hint
from bioflow.configs.main_configs import hint_csv_path
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def get_uniprots_for_hint():
    """
    Recovers UP Gene names maps to UNIPROT nodes containing them.

    :return:
    """
    initial_dict = {}
    for node in DatabaseGraph.get_all('UNIPROT'):
        initial_dict[node._properties['legacyId']] = node.id

    for key in list(initial_dict.keys()):
        initial_dict[key.split('_')[0]] = initial_dict.pop(key)
    return initial_dict


def cross_ref_hint():
    """
    Pulls Hint relationships and connects reached_uniprots_neo4j_id_list in the database

    :return:
    """
    relations_dict = parse_hint(hint_csv_path)
    uniprot_ref_dict = get_uniprots_for_hint()

    processed_nodes = set()
    actual_cross_links = 0
    breakpoints = 300
    size = len(relations_dict)

    log.info('Starting inserting HINT for %s primary nodes' % size)

    for i, (legacyId, linked_legacyIds) in enumerate(relations_dict.items()):

        if i % breakpoints:
            log.info('\t %.2f %%' % (float(i) / float(size) * 100))

        if legacyId in list(uniprot_ref_dict.keys()):
            for linked_legacyId in linked_legacyIds:
                if linked_legacyId in list(uniprot_ref_dict.keys()):
                    actual_cross_links += 1

                    DatabaseGraph.link(uniprot_ref_dict[legacyId], uniprot_ref_dict[linked_legacyId],
                                       'is_interacting', {'source': 'HINT'})


    log.info('HINT Cross-links: %s, HINT processed nodes: %s',
             actual_cross_links, len(processed_nodes))
