"""
Useful thin wrapper of some useful biological knowledge database routine
"""
import os
import pickle
from collections import defaultdict
from csv import reader, writer
from pprint import PrettyPrinter, pprint
from typing import Dict

from bioflow.configs.main_configs import forbidden_neo4j_ids, Dumps
from bioflow.configs.main_configs import reactome_forbidden_nodes, uniprot_forbidden_nodes
from bioflow.utils.general_utils import _is_int
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def look_up_annotation_set(p_load_list, p_type=''):
    """
    Looks up an set of annotations in the database and finds the Ids of nodes containing SWISSPROT
    proteins linked to by annotations

    :param p_load_list: list of payloads
    :param p_type: expected type of payloads
    :return: list of tags for which nodes were not found, list of tags to node property yalues,
    tags to db ids of the first match (preferntially uniprots)
    """
    def db_id_mapping_helper(mapped_db_id_list):
        if mapped_db_id_list:
            return mapped_db_id_list[0][2]
        else:
            return ''

    def transform_annotation_nodes(neo4j_native_nodes):
        retlist = []
        for node in neo4j_native_nodes:
            node_bulbs_id = node.id
            node_legacy_id = node['legacyID']
            node_type = list(node.labels)[0]
            node_display_name = node['displayName']
            payload = (node_type, node_display_name, node_bulbs_id, node_legacy_id)
            if node_type == 'UNIPROT':
                retlist.insert(0, payload)
            else:
                retlist.append(payload)
        return retlist

    load_2_name_list = [(p_load, transform_annotation_nodes(p_nodes))
                        for (p_load, p_nodes) in
                        zip(p_load_list,
                                DatabaseGraph.batch_retrieve_from_annotation_tags(p_load_list,
                                                                                  p_type))]

    db_id_list = [db_id_mapping_helper(value) for key, value in load_2_name_list]
    not_found_list = [key for key, value in load_2_name_list if value == []]
    log.debug('%s IDs out of %s have not been found', len(not_found_list), len(p_load_list))
    log.debug('IDs of missing proteins: %s', not_found_list)
    return not_found_list, load_2_name_list, db_id_list


def _auxilary_annotation_ids_from_csv(source_csv):
    """
    Recovers the set of annotation ids from a csv document

    :param source_csv:
    :return:
    """
    old_ids_parse = []
    new_ids_parse = []

    delimiter = ','

    with open(source_csv, 'rt') as src:
        src_read = src.read()
        # print(src_read)
        if '\t' in src_read:
            delimiter = '\t'
            log.info('tab detected in the %s, switching to tab-separated parse.' % source_csv)

    with open(source_csv, 'rt') as src:
        csv_reader = reader(src, delimiter=delimiter)
        for row in csv_reader:
            if len(row) == 1:
                old_ids_parse = old_ids_parse + row
            if len(row) == 2:
                new_ids_parse.append(row)
            if len(row) > 2:
                raise Exception('more than two values provided in the source file')

    log.debug('debug auxilary_annotation_ids_from_csv: old_ids_parse: %s\n'
             '\tnew_ids_parse:%s' % (old_ids_parse, new_ids_parse))

    if len(old_ids_parse) and len(new_ids_parse):
        log.warning('Both weighted and unweighted hits have been found. Merging and reverting to '
                 'unweighted analysis.')
        old_ids_parse += [_id for _id, _weight in new_ids_parse]

    log.debug('debug: return state old_ids_parse: %s\n'
             '\tnew_ids_parse:%s' % (old_ids_parse, new_ids_parse))

    return old_ids_parse, new_ids_parse


def cast_external_refs_to_internal_ids(analysis_set_csv_location):
    """
    Unwraps the bioflow file specified in the analysis_protein_ids_csv, translates its database
    internal ids for further use

    :param analysis_set_csv_location:
    :return:
    """
    old_live_ids, new_live_ids = _auxilary_annotation_ids_from_csv(analysis_set_csv_location)

    log.debug('debug: \n\tgot old_live_ids (%s)\n\tand new_live_ids (%s)\n\t'
             'from auxilary annotation ids at %s' %
             (old_live_ids, new_live_ids, analysis_set_csv_location))

    if len(old_live_ids) > 0:
        _, _, db_ids_list = look_up_annotation_set(old_live_ids)
        db_ids_list = [_id for _id in db_ids_list if _id != '']
        # # This is not exactly needed and is a more of a log/debug step
        # PrettyPrinter(indent=4, stream=open(Dumps.analysis_set_display_names, 'wt')).pprint(source[1])
        log.debug('old_live_ids branch: mapped ids to %s' % db_ids_list)

        return db_ids_list

    else:
        _, _, db_ids_list = look_up_annotation_set([_id for _id, _weight in new_live_ids])
        log.info('new_live_ids branch: mapped ids to %s' % db_ids_list)
        weighted_db_map = [[_id, float(new_live_ids[i][1])] for i, _id
                           in enumerate(db_ids_list) if _id != '']
        log.info('built the weighted_db_map %s' % weighted_db_map)

        return weighted_db_map

        # writer(open(Dumps.analysis_set_bulbs_ids, 'wt'), delimiter='\n').writerow(weighted_db_map)


def cast_background_set_to_bulbs_id(background_set_csv_location,
                                    analysis_set_csv_location):
    """
    Unwraps the bioflow file specified in the background_internal_ids, translates it to the database
    internal ids for further use

    :param analysis_set_csv_location:
    :param background_set_csv_location:
    :return:
    """
    old_background_ids = []
    if background_set_csv_location:
        old_background_ids, new_background_ids =_auxilary_annotation_ids_from_csv(
            background_set_csv_location)
        old_background_ids += [_id for _id, _weight in new_background_ids]
        old_live_ids, new_live_ids = _auxilary_annotation_ids_from_csv(
            analysis_set_csv_location)
        old_background_ids += old_live_ids
        old_background_ids += [_id for _id, _weight in new_live_ids]
        old_background_ids = list(set(ret for ret in old_background_ids))

    source = look_up_annotation_set(old_background_ids)
    writer(open(Dumps.background_set_bulbs_ids, 'wt'), delimiter='\n').writerow(source[2])


def excluded_nodes_ids_from_names_list():
    """
    Recomputes the list of nodes that contain overloaded terms.

    Without eliminating those terms, they would bring too close together the reactions
    that are normally not, just because of  participation of ultra-abundant elements,
    such as H2O, H+ or ATP

    :param forbidden_entities_list: Dictionary mapping the names of entities to
    their corresponding bulbs classes
    """
    combined_forbidden_set = reactome_forbidden_nodes+uniprot_forbidden_nodes
    forbidden_nodes = DatabaseGraph.mark_forbidden_nodes(combined_forbidden_set)
    forbidden_ids_list = [node.id for node in forbidden_nodes]

    pickle.dump(forbidden_ids_list, open(Dumps.Forbidden_IDs, 'wb'))


def run_diagnostics() -> None:
    """

    :param instructions_list:
    :return:
    """
    DatabaseGraph.node_stats()


def convert_to_internal_ids(base):
    """"
    Converts ids attached to proteins and physical entities to internal db ids, preferentially
    matching to UNIPROTS.

    :param base:
    :return:
    """

    warn_list, results_tuple_list, results_list = look_up_annotation_set(set(base))
    return_dict = {}

    breakpoints = 300
    size = len(results_tuple_list)

    for i, (_key, match_list) in enumerate(results_tuple_list):
        if i % breakpoints == 0:
            # TODO: [loading bar]
            log.info("\t %.2f %%" % (float(i) / float(size)*100))
        if _key not in warn_list:
            for match in match_list:
                if match[0] == 'UNIPROT':
                    return_dict[_key] = match[2]
                else:
                    return_dict[_key] = match_list[0][2]

    log.debug('ID cast converter length: %s', len(return_dict))

    return return_dict


# NOFX: this has to be parametrized and moved to the configs .yaml file for the user to
#  modify the identifiers on which the cross-linking is done.
def cross_link_identifiers():
    """
    Connects all the nodes that share the same uniprot accession numbers. names or gene names
    with an `is_likely_same`

    :return:
    """
    log.info('Cross-linking the identifiers: Uniprot Acnums')
    DatabaseGraph.cross_link_on_xrefs('UNIPROT_Accnum')
    log.info('Cross-linking the identifiers: Uniprot Names')
    DatabaseGraph.cross_link_on_xrefs('UNIPROT_Name')
    log.info('Cross-linking the identifiers: Uniprot Gene Names')
    DatabaseGraph.cross_link_on_xrefs('UNIPROT_GeneName')


def compute_annotation_informativity():
    """
    Computes and writes the number of UNIOPROT nodes each GO term annotates, directly and
    indirectly, computes the informativity of each GO term and finally the total annotation
    information available on the UNIPROT

    :return:
    """
    log.info('Computing the annotation information contents')
    DatabaseGraph.count_go_annotation_cover()


# REFACTOR: [structural diagnostics]: move to structural diagnostics
#   use tabulate for printing the list
#   return the compiled list itself
def pull_up_inf_density():
    """
    Prints out the UIPROT nodes sorted by the amount of information available about them.

    :return:
    """
    name_maps = DatabaseGraph.get_preferential_gene_names()
    print("rank \t informativity \t UNIPROT ID \t gene name")
    for i, node in enumerate(sorted(DatabaseGraph.get_all('UNIPROT'),
                                    key=lambda nde: nde._properties.get('total_information', 0),
                                    reverse=True)):
        print("%4.d \t %.2f \t %s \t %s" % (i+1,
                                            node._properties.get('total_information', 0),
                                            node['legacyID'],
                                            name_maps.get(node['legacyID'], None)))


# REFACTOR: [MAINTAINABILITY] find a more safe and permanent way to do it
on_rtd = os.environ.get('READTHEDOCS') == 'True'
on_unittest = os.environ.get('UNITTESTING') == 'True'

if on_unittest:
    import sys
    import unittests.Mocks.DB_IO_Mocks as SelfMock
    sys.modules[__name__] = SelfMock


def translate_reweight_dict(reweight_dict: Dict) -> Dict:
    """
    Checks a dict assigning desired weight changes in the matrix and if needed translates the
    exterrnal db xrefs to internal db ids

    :param reweight_dict: dict of instructions to re-assign weights
    :return: reweigt_dict translated into internal ids
    """
    discordance_switch = False
    updated_reweight_dict = {}

    for _id_or_tuple, value in reweight_dict.items():
        if type(_id_or_tuple) == tuple:
            if _is_int(_id_or_tuple[0]):
                if discordance_switch:
                    raise Exception("mixed ID types in reweight dict, found an int in %s"
                                    % str(_id_or_tuple))
                return reweight_dict  # internal DB ids were supplied
            else:
                discordance_switch = True
                id_lookup = convert_to_internal_ids(_id_or_tuple)
                updated_reweight_dict[(id_lookup[_id_or_tuple[0]],
                                       id_lookup[_id_or_tuple[1]])] = float(value)

        else:
            if _is_int(_id_or_tuple):
                if discordance_switch:
                    raise Exception("mixed ID types in reweight dict, found an int in %s"
                                    % str(_id_or_tuple))
                return reweight_dict
            else:
                discordance_switch = True
                log.debug('debug of c: %s' % _id_or_tuple)
                id_lookup = convert_to_internal_ids([_id_or_tuple])
                log.debug('id lookup parse: %s' % str(id_lookup))
                updated_reweight_dict[(id_lookup[_id_or_tuple])] = float(value)

    return updated_reweight_dict

if __name__ == "__main__":
    # erase_custom_fields()
    # excluded_nodes_ids_from_names_list(forbidden_verification_list)
    # print cast_external_refs_to_internal_ids()
    # print look_up_annotation_set('CTR86')
    # print look_up_annotation_set('ENSG00000131981', 'UNIPROT_Ensembl')
    # cast_external_refs_to_internal_ids()
    # cast_background_set_to_bulbs_id()

    # _, resdict, reslist = look_up_annotation_set(['RNF14'])
    # pprint(resdict)
    # print reslist
    # _, resdict, reslist = look_up_annotation_set(['RNF14'])
    # pprint(resdict)
    # print reslist

    # compute_annotation_informativity()
    # pull_up_inf_density()

    # run_diagnostics(to_deprecate_full_list)
    # memoize_bulbs_type(to_deprecate_neo4j_names_dict['UNIPROT'][0])
    # cast_external_refs_to_internal_ids()
    # cast_background_set_to_bulbs_id(background_set_csv_location=None)
    # Akshay_p53_go_set = ["0030330", "0000019", "0000002", "0006977"]
    # go_bulbs_ids = lookup("GO Term", Akshay_p53_go_set)
    # association_dict = pull_associations(go_bulbs_ids, ["is_go_annotation"])
    # proteins = [item for key, val in association_dict.iteritems() for item in val]
    # writer(open(Dumps.analysis_set_bulbs_ids, 'w'), delimiter='\n').writerow(proteins)
    # print proteins

    bulbs_id = look_up_annotation_set(['FAA4'])
    print(bulbs_id)
