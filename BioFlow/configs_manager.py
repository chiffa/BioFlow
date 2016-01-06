"""
:created: Dec 15, 2013
:author: Andrei Kucharavy
Performs IO from the setup .ini files and casts into relevant Python Dictionaries
"""
from os.path import join, abspath, expanduser
import os
from string import lower
from BioFlow.utils.general_utils.dict_like_configs_parser import ini_configs2dict, dict2init_configs
from BioFlow.utils.general_utils.high_level_os_io import mkdir_recursive
from BioFlow.utils.general_utils.internet_io import url_to_local
from collections import defaultdict


configs_rootdir = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        '../configs/'))

shortnames = ['servers', 'options', 'sources', 'predictions']
configsfiles = dict([(name, join(configs_rootdir, name + '.ini'))
                     for name in shortnames])


class StructureGenerator(object):

    _online_DBs = {
        'REACTOME': [r'http://www.reactome.org/download/current/biopax.zip'],
        'UNIPROT': [r'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'],
        'HINT': [
            r'http://hint.yulab.org/CervBinaryHQ.txt',
            r'http://hint.yulab.org/HumanBinaryHQ.txt',
            r'http://hint.yulab.org/MouseBinaryHQ.txt'],
        'GO': [r'http://purl.obolibrary.org/obo/go/go-basic.obo'],
        'BIOGRID': [r'http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.3.122/BIOGRID-ORGANISM-3.3.122.tab2.zip'],
        'ABOUNDANCES': []}

    # paths to be appended to the user-provided installation directory
    _local_file_tree = {
        'INTERNAL': '.',
        'REACTOME': 'Reactome',
        'UNIPROT': 'Uniprot/uniprot_sprot.dat',
        'HINT': 'HiNT',
        'GO': 'GO/go.obo',
        'BIOGRID': 'BioGRID',
        'SIDER': 'SIDER2/meddra_adverse_effects.tsv',
        'ABOUNDANCES': 'Protein_aboundances',
        'CHROMOSOMES': 'Chr_mappings'}

    # default configuration elements for yeast protein analysis
    _S_Cerevisae = {
        'shortname': 'sCerevisae',
        'tax_id': '559292',
        'Reactome_name': 'Saccharomyces_cerevisiae.owl',
        'Biogrid_name': 'Saccharomyces_cerevisiae',
        'HINT_name': 'CervBinaryHQ.txt'}

    # default configuration elements for human proteins analysis
    _Human = {
        'shortname': 'human',
        'tax_id': '9606',
        'Reactome_name': 'Homo_sapiens.owl',
        'Biogrid_name': 'Homo_Sapiens',
        'HINT_name': 'HumanBinaryHQ.txt'}

    # default configuration elements for mice proteins analysis
    _Mice = {
        'shortname': 'mouse',
        'tax_id': '10090',
        'Reactome_name': 'Mus_musculus.owl',
        'Biogrid_name': 'Mus_musculus',
        'HINT_name': 'HumanBinaryHQ.txt'}

    reforgs = ['mouse', 'human', 'yeast']

    @classmethod
    def _generate_template(cls, payload_dict, expanded=False):
        """
        Generates a template dictionary that would be converted to a sources.ini config file
            to a file that would be
        :param payload_dict: Dict containing tax_id, Reactome and biogrid filenames. if expanded
        use is anticipated, needs chromosome namepattern
        :param expanded: if true, will try to add additional options into the configs file
        :return:
        """
        template_dict = {
            'INTERNAL': {'mongoprefix': '_' + payload_dict['shortname'],
                         'mongosuffix': '_v_1',
                         'dumpprefix': '/' + payload_dict['shortname'],
                         'load': 'NotAFile.txt'},
            'REACTOME':
                {'load': payload_dict['Reactome_name']},
            'UNIPROT':
                {'tax_ids': payload_dict['tax_id']},
            'HINT':
                {'load': payload_dict['HINT_name']},
            'BIOGRID':
                {'load': payload_dict['Biogrid_name']},
            'GO':
                {}
        }

        if expanded:
            additional_options = {
                'CHROMOSOMES': {'load': payload_dict['name_pattern'],
                                'namepattern': payload_dict['name_pattern']},
                'ABOUNDANCES': {'load': payload_dict['tax_id']},
            }
        else:
            additional_options = {
                'CHROMOSOMES': {'load': '-1', 'namepattern': '-1'},
                'ABOUNDANCES': {'load': '-1'},
            }
        template_dict.update(additional_options)

        return template_dict

    @classmethod
    def _add_location_to_template(cls, template_dict):
        """
        Adds a location to the template dictionary used to generate a configuration file
        :param template_dict:
        :return:
        """
        master_location = ini_configs2dict(configsfiles['servers'])[
            'PRODUCTION']['base_folder']
        for key, value in template_dict.iteritems():
            value['location'] = join(
                master_location, cls._local_file_tree[key])
        return template_dict

    @classmethod
    def build_source_config(cls, pl_type):
        """
        Writes a source file based on the string organism argument

        :param pl_type: string falling into one of the three categories. if not valid, an custom
        exception is raised
        :return:
        """
        if pl_type not in cls.reforgs:
            raise Exception(
                'Unsupported organism, %s not in %s.' %
                (pl_type, cls.reforgs) + ' Please modify the sources.ini manually')
        else:
            write_path = join(configs_rootdir, 'sources.ini')
            cfdict = {}
            if pl_type == 'mouse':
                cfdict = cls._generate_template(cls._Mice)
            if pl_type == 'human':
                cfdict = cls._generate_template(cls._Human)
            if pl_type == 'yeast':
                cfdict = cls._generate_template(cls._S_Cerevisae)
            cfdict = cls._add_location_to_template(cfdict)
            print cfdict
            dict2init_configs(write_path, cfdict)

    @classmethod
    def pull_online_dbs(cls):
        """
        Pulls the databases mentionned online to the direction in the

        :return:
        """
        # read the sources.ini file
        write_dirs = ini_configs2dict(configsfiles['sources'])
        write_dirs = dict([(key, directory['location']) for key, directory in
                           write_dirs.iteritems() if key not in ['INTERNAL', 'CHROMOSOMES']])
        for DB_type, location in write_dirs.iteritems():
            # recursively create directories if the directories don't exist
            for sublocation in cls._online_DBs[DB_type]:
                print 'loading %s database from %s to %s' % (lower(DB_type), sublocation, location)
                mkdir_recursive(location)
                url_to_local(sublocation, location)


def graceful_fail():
    raise Exception('Argument not contained in the initial config file')


def parse_configs():
    """ parses all the relevant configs and inserts graceful failure to all of them """
    pre_srces = ini_configs2dict(configsfiles['sources'])
    srces = defaultdict(graceful_fail)
    srces.update(pre_srces)

    return ini_configs2dict(configsfiles['servers']), \
        ini_configs2dict(configsfiles['options']),\
        srces,\
        ini_configs2dict(configsfiles['predictions'])


def conf_file_path_flattener(raw_configs_dict):
    """
    Rips out the 'location' argument and parses it into a filename from the parse of the
    configurations tree and flattens the resulting dictionary. It is used mainly for the
    accelerated file location recovery

    :param raw_configs_dict: Initial file location and parameter locations: `
    {'HINT':{'load': ld, 'location': loc}...}`
    :return: `{'HINT': }`
    """
    paths_dict = {}
    for ext_DB_type, ext_DB_args in raw_configs_dict.iteritems():
        current_path = ext_DB_args['location']
        if os.path.isdir(current_path):
            for fle in [f for f in os.listdir(current_path) if
                        os.path.exists(os.path.join(current_path, f))]:
                if ext_DB_args['load'] in fle:
                    current_path = os.path.join(current_path, fle)
        paths_dict[ext_DB_type] = current_path
    final_paths_dict = defaultdict(graceful_fail)
    final_paths_dict.update(paths_dict)
    return final_paths_dict


def edit_confile(conf_shortname, section, parameter, newvalue):
    """

    :param conf_shortname:
    :param section:
    :param parameter:
    :param newvalue:
    :return:
    """
    tmp_confdict = ini_configs2dict(configsfiles[conf_shortname])
    tmp_confdict[section][parameter] = newvalue
    dict2init_configs(configsfiles[conf_shortname], tmp_confdict)


def set_folders(file_directory, neo4jserver='http://localhost:7474',
                mongoserver='mongodb://localhost:27017/'):
    if file_directory[0] == r'~':
        file_directory = expanduser(file_directory)
    edit_confile(
        'servers',
        'PRODUCTION',
        'base_folder',
        abspath(file_directory))
    edit_confile('servers', 'PRODUCTION', 'server_neo4j', neo4jserver)
    edit_confile('servers', 'PRODUCTION', 'mongodb_server', mongoserver)
    StructureGenerator.build_source_config('yeast')


if __name__ == "__main__":
    set_folders('/home/ank/datastore')
    StructureGenerator.build_source_config('yeast')
    # pp = PrettyPrinter(indent=4)
    # pp.pprint(parse_configs())
    # StructureGenerator.pull_online_dbs()
