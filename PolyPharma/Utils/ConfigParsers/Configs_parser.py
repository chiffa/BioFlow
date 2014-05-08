'''
Created on Dec 15, 2013
@author: ank
performs IO
'''

from __future__ import print_function
import ConfigParser
from Configs_common import configsfiles
from pprint import PrettyPrinter
import os

def parse_configs():
    ''' parses all the relevant configs '''

    def improved_read(path,):
        '''
        Parses a config file on given path, in case of failure raises an IOError, 
        '''
        cfg=ConfigParser.SafeConfigParser()
        rfs=cfg.read(path)
        if rfs==[]:
            raise IOError('cannot load '+path)
        MainDict={}
        for section in cfg.sections():
            section_name="_".join([elt for elt in section.split('_')])
            MainDict[section_name]={}
            for option in cfg.options(section):
                MainDict[section_name][option]=cfg.get(section, option)
        return MainDict
    
    return improved_read(configsfiles[0]),\
           improved_read(configsfiles[1]),\
           improved_read(configsfiles[2]),\
           improved_read(configsfiles[3]),

def sourcefile_compilator(Sources_dict):
    finalPathdict = {}
    for ext_DB_type, ext_DB_args in Sources_dict.iteritems():
        currentPath = ext_DB_args['location']
        if os.path.isdir(currentPath):
            for fle in [f for f in os.listdir(currentPath) if os.path.exists(os.path.join(currentPath,f))]:
                if ext_DB_args['load'] in fle:
                    currentPath = os.path.join(currentPath, fle)
        finalPathdict[ext_DB_type] = currentPath
    return finalPathdict


if __name__ == "__main__":
    pp=PrettyPrinter(indent=4)
    pp.pprint(parse_configs())