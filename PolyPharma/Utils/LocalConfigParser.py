"""
:created: Dec 15, 2013
:author: Andrei Kucharavy
Performs IO from the setup .ini files and casts into relevant Python Dictionaries
"""
from pprint import PrettyPrinter
from os.path import join
import os
from PolyPharma.Utils.GeneralUtils.SanerConfigsParser import ini_configs2dict, dict2init_configs
from collections import defaultdict

configs_rootdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../configs/'))

shortnames = ['servers', 'options', 'sources', 'predictions']
configsfiles = dict([(name, join(configs_rootdir, name+'.ini')) for name in shortnames ])


class StructureGenerator(object):

    online_file_tree = { }

    # paths to be appended to the user-provided installation directory
    local_file_tree = { 'META': '-1',
                        'REACTOME' : 'Reactome',
                        'UNIPROT' : 'Uniprot/uniprot_sprot.dat',
                        'HINT' : 'HiNT',
                        'GO' : 'GO',
                        'BIOGIRD' : 'BioGRID',
                        'SIDER': 'SIDER2/meddra_adverse_effects.tsv',
                        'ABOUNDANCES': 'Protein_aboundances',
                        'CHROMOSOMES': 'Chr_mappings'}

    # default configuration elements for yeast protein analysis
    S_Cerevisae = {'shortname': 'sCerevisae',
                   'tax_id': '559292',
                   'Reactome_name': 'Saccharomyces cerevisiae.owl',
                   'Biogrid_name': 'Saccharomyces_cerevisae.tsv'}

    # default configuration elements for human proteins analysis
    Human = {      'shortname': 'human',
                   'tax_id': '9606',
                   'Reactome_name': 'Homo Sapiens.owl',
                   'Biogrid_name': 'Homo_Sapiens.tsv'}

    # default configuration elements for mice proteins analysis
    Mice = {       'shortname': 'mouse',
                   'tax_id': '10090',
                   'Reactome_name': 'Mus Musculus.owl',
                   'Biogrid_name': 'Mus_musculus.tsv'}

    @classmethod
    def generate_template(cls, payload_dict, expanded=False):
        """
        Generates a template dictionary that would be converted to a sources.ini config file
            to a file that would be
        :param payload_dict: Dict containing tax_id, Reactome and biogrid filenames. if expanded use is anticipated, needs chromosome namepattern
        :param expanded: if true, will try to add additional options into the configs file
        :return:
        """
        template_dict = {'META':{'mongoprefix': '_'+payload_dict['shortname'],
                                 'mongosuffix': '_v_1',
                                 'dumpprefix': '/'+payload_dict['shortname']},
                                'REACTOME' :
                                    {'load': payload_dict['Reactome_name']},
                                'UNIPROT':
                                    {'tax_ids': payload_dict['tax_id']},
                                'HINT':
                                    {'load': payload_dict['tax_id']},
                                'BIOGIRD':
                                    {'load': payload_dict['Biogrid_name']},
                               }

        if expanded:
            additional_options =  {
                                'CHROMOSOMES':
                                    {'load': payload_dict['name_pattern'],
                                     'namepattern': payload_dict['name_pattern']},
                                'ABOUNDANCES':
                                    {'load': payload_dict['tax_id']},
                                }
        else:
            additional_options = {
                                'CHROMOSOMES':
                                    {'load': '-1',
                                     'namepattern': '-1'},
                                'ABOUNDANCES':
                                    {'load': '-1'},
                                }
        template_dict.update(additional_options)

        return template_dict

    @classmethod
    def add_location_to_template(cls, template_dict):
        """
        Adds a location to the template dictionary used to generate a configuration file
        :param template_dict:
        :return:
        """
        master_location = ini_configs2dict(configsfiles['options'])['MASTER']['base_folder']
        for key, value in template_dict.iteritems():
            value['location'] = join(master_location, cls.local_file_tree[key])
        return template_dict

    @classmethod
    def build_source_config(cls, pl_type):
        """
        Writes a source file based on the string organism argument

        :param pl_type: string falling into one of the three categories. if not valid, an costum exception is raised
        :return:
        """
        if pl_type not in ['mouse', 'human', 'yeast']:
            raise Exception('Unsupported organism, not in %s. Please modify the sources.ini manually' % pl_type)
        else:
            write_path = join(configs_rootdir, 'test_sources.ini')
            cfdict = {}
            if pl_type == 'mouse':
                cfdict = cls.generate_template(cls.Mice)
            if pl_type == 'human':
                cfdict = cls.generate_template(cls.Human)
            if pl_type == 'yeast':
                cfdict = cls.generate_template(cls.S_Cerevisae)
            cfdict = cls.add_location_to_template(cfdict)
            print cfdict
            dict2init_configs(write_path, cfdict)



def graceful_fail():
    raise Exception('Argument not contained in the initial config file')


def parse_configs():
    """ parses all the relevant configs and inserts graceful failure to all of them """
    pre_srces = ini_configs2dict(configsfiles['sources'])
    srces = defaultdict(graceful_fail )
    srces.update(pre_srces)

    return ini_configs2dict(configsfiles['servers']),\
           ini_configs2dict(configsfiles['options']),\
           srces,\
           ini_configs2dict(configsfiles['predictions']),


def conf_file_path_flattener(raw_configs_dict):
    """
    Rips out the 'location' argument and parses it into a filename from the parse of the configurations tree and flattens the resulting dictionary
    Is used mainly for the accelerated file location recovery

    :param raw_configs_dict: Initial file location and parameter locations: `{'HINT':{'load': ld, 'location': loc}...}`
    :return: `{'HINT': }`
    """
    paths_dict = {}
    for ext_DB_type, ext_DB_args in raw_configs_dict.iteritems():
        currentPath = ext_DB_args['location']
        if os.path.isdir(currentPath):
            for fle in [f for f in os.listdir(currentPath) if os.path.exists(os.path.join(currentPath, f))]:
                if ext_DB_args['load'] in fle:
                    currentPath = os.path.join(currentPath, fle)
        paths_dict[ext_DB_type] = currentPath
    final_paths_dict = defaultdict(graceful_fail)
    final_paths_dict.update(paths_dict)
    return final_paths_dict


if __name__ == "__main__":
    StructureGenerator.build_source_config('yeast')
    # pp = PrettyPrinter(indent=4)
    # pp.pprint(parse_configs())
    pass