"""
Advanced logic functions required for configurations to load
"""
import os
from typing import Tuple
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def compute_full_paths(sources_parse_dict: dict,
                       online_dbs_parse_dict: dict,
                       base_folder: dict) -> dict:
    """
    Computes all of the base files locations for the set or sources, based on the locations
    specified in the online_dbs.inin configuration file

    :param sources_parse_dict:
    :param online_dbs_parse_dict:
    :param base_folder:
    :return:
    """
    paths_dict = {}

    for source_name, source_contents in list(sources_parse_dict.items()):

        if source_name in ['INTERNAL', 'META']:
            continue

        source_is_inactive = bool(online_dbs_parse_dict[source_name].get('inactive', False))

        if not online_dbs_parse_dict[source_name] or source_is_inactive:
            log.exception('Attempted to load an inactive or inexistent source type %s', source_name)
            continue

        pre_location = online_dbs_parse_dict[source_name]['local'].replace('$DB_HOME$', base_folder)

        if 'file' in list(source_contents.keys()):
            paths_dict[source_name] = os.path.join(pre_location, source_contents['file'])

        elif 'name_pattern' in list(source_contents.keys()):
            target = ''

            if not os.path.exists(pre_location):
                target = 'Invalid for now'
                log.warning('name_pattern %s was not matched in %s. Pull the online dbs first',
                            source_contents['name_pattern'], pre_location)

            else:
                for file_name in os.listdir(pre_location):
                    # Now raises an error if the folder has not been initialized
                    if source_contents['name_pattern'] in file_name:
                        target = file_name
                        break

            if target:
                paths_dict[source_name] = os.path.join(pre_location, target)

            else:
                log.exception('cannot find name pattern %s in folder %s',
                              source_contents['name_pattern'],
                              pre_location)

        else:
            paths_dict[source_name] = pre_location

    return paths_dict


if __name__ == "__main__":
    pass
    # pull_online_dbs()
    # pp = PrettyPrinter(indent=4)
    # pp.pprint(parse_configs())
    # pull_online_dbs()
