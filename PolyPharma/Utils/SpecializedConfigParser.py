"""
:created: Dec 15, 2013
:@author: Andrei Kucharavy
performs IO from the setup .ini files and casts into relevant Python Dictionaries
"""
from pprint import PrettyPrinter
from os.path import join
import os
from PolyPharma.Utils.GeneralUtils.SanerConfigsParser import ini_configs2dict
from collections import defaultdict

configs_rootdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../configs/'))

shortnames = ['servers', 'options', 'sources', 'predictions']
configsfiles = dict([(name, join(configs_rootdir, name+'.ini')) for name in shortnames ])


def graceful_fail():
    raise Exception('custom_exception')


def parse_configs():
    """ parses all the relevant configs """
    pre_srces = ini_configs2dict(configsfiles['sources'])
    srces = defaultdict(graceful_fail )
    srces.update(pre_srces)

    return ini_configs2dict(configsfiles['servers']),\
           ini_configs2dict(configsfiles['options']),\
           srces,\
           ini_configs2dict(configsfiles['predictions']),


def conf_file_path_flattener(raw_configs_dict):
    """
    Rips out the 'location' argument from the parse of the configurations tree and flattens the resulting dictionary
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
    pp = PrettyPrinter(indent=4)
    pp.pprint(parse_configs())