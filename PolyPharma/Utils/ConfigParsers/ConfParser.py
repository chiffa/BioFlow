'''
Created on Dec 15, 2013
@author: ank
performs configfiles reading and writing
'''

# Python 2.7 specific imports here
from __future__ import print_function
import ConfigParser
import os
from pprint import PrettyPrinter

rootdir=os.path.abspath(os.path.join(os.path.dirname(__file__),'../configs/'))
shortnames=['servers','options']
configsfiles=[rootdir+'/'+name+'.ini' for name in shortnames ]

def generate_configs():
    '''
    Generates the configFiles. Executed only on call of module as if it was main
    '''
    
    def generate_servers_config():
        '''
        Generates configuration files for the different execution types
        '''
        config = ConfigParser.SafeConfigParser()
        config.set('DEFAULT', "servers_are_local", "True")
        config.set('DEFAULT', "server_neo4j","http://localhost:7474")
        config.add_section('TEST')
        config.set('TEST', "local_neo4j","/Users/ank/Programming/DBs/neo4j")
        config.add_section('PRODUCTION')
        config.set('PRODUCTION', "local_neo4j","/usr/local/Cellar/neo4j/1.9.5/libexec/data")
        with open(configsfiles[0],'w') as configfile:
            config.write(configfile)
    
    def generate_options_config():
        '''
        Generates configurations files for different runtime options
        '''
        config = ConfigParser.SafeConfigParser()
        config.add_section('JVM')
        config.set('JVM',"NEO4J_PYTHON_JVMARGS", "-Xms128M -Xmx512M")
        config.set('JVM',"JAVA_HOME", "/usr/lib/jvm/java")
        config.add_section('LOGGER')
        config.set('LOGGER',"level","DEBUG")
        config.set('LOGGER',"output_file","/Users/ank/Programming/logs/production_neo4j_logger.txt")
        
        with open(configsfiles[1],'w') as configfile:
            config.write(configfile)
    
    generate_options_config()
    generate_servers_config()

def parse_configs():
    ''' parses all the relevant configs '''

    def improved_read(path):
        '''
        Parses a config file on given path, in case of failure raises an IOError, 
        '''
        cfg=ConfigParser.SafeConfigParser()
        rfs=cfg.read(path)
        if rfs==[]:
            raise IOError('cannot load '+path)
        MainDict={}
        for section in cfg.sections():
            section_name="_".join([elt.lower() for elt in section.split('_')])
            MainDict[section_name]={}
            for option in cfg.options(section):
                MainDict[section_name][option]=cfg.get(section, option).lower()
        return MainDict
    
    return improved_read(configsfiles[0]), improved_read(configsfiles[1])


if __name__ == "__main__":
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)
    generate_configs()
    pp=PrettyPrinter(indent=4)
    pp.pprint(parse_configs())