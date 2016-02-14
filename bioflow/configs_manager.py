"""
:created: Dec 15, 2013
:author: Andrei Kucharavy
Performs IO from the setup .ini files and casts into relevant Python Dictionaries
"""
from os.path import join, abspath, expanduser
import os
from shutil import copy as copy_file
from string import lower
from bioflow.utils.general_utils.dict_like_configs_parser import ini_configs2dict, dict2init_configs
from bioflow.utils.general_utils.high_level_os_io import mkdir_recursive
from bioflow.utils.general_utils.internet_io import url_to_local
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)

configs_rootdir = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        'configs/'))

conf_files_shortnames = ['servers', 'options', 'sources', 'predictions', 'online_dbs']
conf_files_locations = dict([(name, join(configs_rootdir, name + '.ini'))
                             for name in conf_files_shortnames])

ref_orgs_shortnames = [
    shortname.strip() for shortname in
    ini_configs2dict(conf_files_locations['options'])['ORGANISMS']['allowed'].split(',')]
ref_orgs_locations = dict([(name, join(configs_rootdir, 'reference', name + '.ini'))
                           for name in ref_orgs_shortnames])


def build_source_config(organism_shortname):
    """
    Writes a bioflow file based on the string organism argument

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


def pull_online_dbs():
    """
    Pulls the databases mentionned online to the direction in the

    :return:
    """
    write_dirs = ini_configs2dict(conf_files_locations['online_dbs'])
    base_folder = ini_configs2dict(conf_files_locations['servers'])['PRODUCTION']['base_folder']

    for DB_type, location_dict in write_dirs.iteritems():

        if 'inactive' in location_dict.keys() and location_dict['inactive'] == 'True':
            continue

        local = location_dict['local'].replace('$DB_HOME$', base_folder)
        onlines = [location.strip() for location in location_dict['online'].split(',')]

        for online in onlines:
            log.info('loading %s database from %s to %s',
                     lower(DB_type), location_dict['online'], local)
            mkdir_recursive(local)
            url_to_local(online, local)


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
    for source_name, source_contents in sources_parse_dict.items():

        if source_name == 'INTERNAL':
            continue

        if not online_dbs_parse_dict[source_name] or \
                ('inactive' in online_dbs_parse_dict[source_name].keys() and
                 online_dbs_parse_dict[source_name]['inactive'] == 'True'):
            log.exception('attempt to load an inactive source type %s', source_name)
            continue

        pre_location = online_dbs_parse_dict[source_name]['local'].replace('$DB_HOME$', base_folder)

        if 'file' in source_contents.keys():
            paths_dict[source_name] = join(pre_location, source_contents['file'])

        elif 'name_pattern' in source_contents.keys():
            target = ''

            for file_name in os.listdir(pre_location):
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
                neo4jserver='http://localhost:7474',
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
    build_source_config('yeast')


if __name__ == "__main__":
    set_folders('/home/andrei/support')
    build_source_config('human')
    # pp = PrettyPrinter(indent=4)
    # pp.pprint(parse_configs())
    # pull_online_dbs()
