__author__ = 'ank'

import ConfigParser
from PolyPharma.Utils.ConfigParsers.Configs_common import rootdir,configsfiles,shortnames
import os

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

    if not os.path.exists(rootdir):
        os.makedirs(rootdir)
    generate_options_config()
    generate_servers_config()

if __name__=="__main__":
    generate_configs()