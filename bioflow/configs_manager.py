"""
:created: Dec 15, 2013
:author: Andrei Kucharavy
Performs IO from the setup .ini files and casts into relevant Python Dictionaries
"""
from os.path import join, abspath, expanduser
import os
from shutil import copy as copy_file
from bioflow.utils.general_utils.dict_like_configs_parser import ini_configs2dict, dict2init_configs
from bioflow.utils.general_utils.high_level_os_io import mkdir_recursive
from bioflow.utils.general_utils.internet_io import url_to_local, marbach_post_proc
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)

configs_rootdir = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        'configs/'))

conf_files_shortnames = ['servers', 'sources', 'online_dbs']
conf_files_locations = dict([(name, join(configs_rootdir, name + '.ini'))
                             for name in conf_files_shortnames])

ref_orgs_shortnames = ['mouse', 'human', 'yeast']
ref_orgs_locations = dict([(name, join(configs_rootdir, 'reference', name + '.ini'))
                           for name in ref_orgs_shortnames])


def build_source_config(organism_shortname):
    """
    Copies BioFlow reference config file for an organism to the active config file

    :param organism_shortname: string falling into one of the three categories. if not valid, an custom
    exception is raised
    :return:
    """
    if organism_shortname not in ref_orgs_shortnames:
        raise Exception('Unsupported organism, %s not in %s. %s',
                        (organism_shortname, ref_orgs_shortnames,
                         'Please modify the sources.ini manually'))

    else:
        write_path = join(configs_rootdir, 'sources.ini')
        copy_file(ref_orgs_locations[organism_shortname], write_path)
        log.info('Active organism set to %s' % organism_shortname)


def pull_online_dbs():
    """
    Pulls the databases mentionned online to the direction in the

    :return:
    """
    write_dirs = ini_configs2dict(conf_files_locations['online_dbs'])
    base_folder = ini_configs2dict(conf_files_locations['servers'])['PRODUCTION']['base_folder']

    for DB_type, location_dict in write_dirs.items():

        print(DB_type)

        if 'inactive' in list(location_dict.keys()) and location_dict['inactive'] == 'True':
            continue

        local = location_dict['local'].replace('$DB_HOME$', base_folder)
        onlines = [location.strip() for location in location_dict['online'].split(',')]

        if 'rename' in list(location_dict.keys()):
            renames = [location.strip() for location in location_dict['rename'].split(',')]

            for online, rename in zip(onlines, renames):
                log.info('loading %s database from %s to %s',
                         DB_type.lower(), location_dict['online'], local)
                mkdir_recursive(local)
                url_to_local(online, local, rename=rename)

        else:
            for online in onlines:
                log.info('loading %s database from %s to %s',
                         DB_type.lower(), location_dict['online'], local)
                mkdir_recursive(local)
                url_to_local(online, local)

        if DB_type == 'MARBACH':
            log.info('cleaning up local marbach databse')
            marbach_post_proc(local)


def parse_config(configfile_shortname):
    """
    Parses the config file given config short name

    :param configfile_shortname:
    :return:
    """
    return ini_configs2dict(conf_files_locations[configfile_shortname])


def compute_full_paths(sources_parse_dict, online_dbs_parse_dict, servers_parse_dict):
    """
    Computes all of the base files locations for the set or sources, based on the locations
    specified in the online_dbs.inin configuration file

    :param sources_parse_dict:
    :param online_dbs_parse_dict:
    :param servers_parse_dict:
    :return:
    """
    base_folder = servers_parse_dict['base_folder']
    paths_dict = {}

    for source_name, source_contents in list(sources_parse_dict.items()):

        if source_name in ['INTERNAL', 'META']:
            continue

        if not online_dbs_parse_dict[source_name] or \
                ('inactive' in list(online_dbs_parse_dict[source_name].keys()) and
                 online_dbs_parse_dict[source_name]['inactive'] == 'True'):
            log.exception('attempt to load an inactive source type %s', source_name)
            continue

        pre_location = online_dbs_parse_dict[source_name]['local'].replace('$DB_HOME$', base_folder)

        if 'file' in list(source_contents.keys()):
            paths_dict[source_name] = join(pre_location, source_contents['file'])

        elif 'name_pattern' in list(source_contents.keys()):
            target = ''

            if not os.path.exists(pre_location):
                target = 'Invalid for now'
                log.warning('name_pattern %s was not matched in %s. Pull the online dbs first',
                            source_contents['name_pattern'], pre_location)

            else:
                for file_name in os.listdir(pre_location):
                    # TODO: now raises an error if the folder has not been initialized
                    if source_contents['name_pattern'] in file_name:
                        target = file_name
                        break

            if target:
                paths_dict[source_name] = join(pre_location, target)

            else:
                log.exception('cannot find name pattern %s in folder %s',
                              source_contents['name_pattern'],
                              pre_location)

        else:
            paths_dict[source_name] = pre_location

    return paths_dict


def edit_config_file(conf_shortname, section, parameter, new_value):
    """

    :param conf_shortname:
    :param section:
    :param parameter:
    :param new_value:
    :return:
    """
    tmp_config_dict = ini_configs2dict(conf_files_locations[conf_shortname])
    tmp_config_dict[section][parameter] = new_value
    dict2init_configs(conf_files_locations[conf_shortname], tmp_config_dict)


def set_folders(file_directory,
                neo4jserver='bolt://localhost:7687',
                mongoserver='mongodb://localhost:27017/'):
    if file_directory[0] == r'~':
        file_directory = expanduser(file_directory)
    edit_config_file(
        'servers',
        'PRODUCTION',
        'base_folder',
        abspath(file_directory))
    edit_config_file('servers', 'PRODUCTION', 'server_neo4j', neo4jserver)
    edit_config_file('servers', 'PRODUCTION', 'mongodb_server', mongoserver)
    build_source_config('human')
    log.info('Download folders and servers configs are set')


if __name__ == "__main__":
    set_folders('/home/andrei/sources')
    build_source_config('human')
    pull_online_dbs()
    # pp = PrettyPrinter(indent=4)
    # pp.pprint(parse_configs())
    # pull_online_dbs()
