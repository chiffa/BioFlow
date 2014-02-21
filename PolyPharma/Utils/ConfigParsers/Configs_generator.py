__author__ = 'ank'

import ConfigParser
from PolyPharma.Utils.ConfigParsers.Configs_common import rootdir, configsfiles, \
                                            db_root, source_ext_db_roots, source_ext_prediction_roots
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
        config.set('DEFAULT', "server_neo4j", "http://localhost:7474")
        config.set('DEFAULT', "local_sqlite", db_root+"Data_store/sqlite")
        config.set('DEFAULT', "mongoDB_server","mongodb://localhost:27017/")
        config.add_section('TEST')
        config.set('TEST', "local_neo4j", db_root+"Data_store/neo4j")
        config.add_section('PRODUCTION')
        config.set('PRODUCTION', "local_neo4j", db_root+"neo4j-human/data")

        with open(configsfiles[0], 'w') as configfile:
            config.write(configfile)

        return config

    def generate_options_config():
        '''
        Generates configurations files for different runtime options
        '''
        config = ConfigParser.SafeConfigParser()
        config.add_section('JVM')
        config.set('JVM', "NEO4J_PYTHON_JVMARGS", "-Xms128M -Xmx512M")
        config.set('JVM', "JAVA_HOME", "/usr/lib/jvm/java")
        config.add_section('LOGGER')
        config.set('LOGGER', "level", "DEBUG")
        config.set('LOGGER', "output_file", db_root+"logs/production_neo4j_logger.txt")

        with open(configsfiles[1], 'w') as configfile:
            config.write(configfile)
        return config

    def generate_source_database_config():
        '''
        Generates configurations for reading from the external databases
        '''
        config = ConfigParser.SafeConfigParser()
        config.add_section('REACTOME')
        config.set('REACTOME', 'location', source_ext_db_roots+'Reactome')
        config.set('REACTOME', 'load', 'Homo sapiens.owl')
        config.add_section('UNIPROT')
        config.set('UNIPROT', 'location', source_ext_db_roots+'Uniprot/uniprot_sprot.dat')
        config.set('UNIPROT', 'Tax_IDs', '9606, ')
        config.add_section('HINT')
        config.set('HINT', 'location', source_ext_db_roots+'HiNT')
        config.set('HINT', 'load', 'HumanBinaryHQ.txt')
        config.add_section('GO')
        config.set('GO', 'location', source_ext_db_roots+'GO/go.obo')
        config.add_section('ABOUNDANCES')
        config.set('ABOUNDANCES', 'location', source_ext_db_roots+'Protein_aboundances')
        config.set('ABOUNDANCES', 'load', '9606')
        config.add_section('SIDER')
        config.set('SIDER', 'location', source_ext_db_roots+'SIDER2/meddra_adverse_effects.tsv')

        with open(configsfiles[2], 'w') as configfile:
            config.write(configfile)
        return config

    def generate_source_predictions_config():
        '''
        Generates configuration files for reading from external prediction datasets
        '''
        config = ConfigParser.SafeConfigParser()
        config.add_section('OVERINGTON')
        config.set('OVERINGTON', 'location', source_ext_prediction_roots+'Overington_raw')
        config.add_section('NEFLANAVIR')
        config.set('NEFLANAVIR', 'location', source_ext_prediction_roots+'NeflanavirSource.csv')
        with open(configsfiles[3], 'w') as configfile:
            config.write(configfile)
        return config

    if not os.path.exists(rootdir):
        os.makedirs(rootdir)
    generate_options_config()
    generate_servers_config()
    generate_source_database_config()
    generate_source_predictions_config()

if __name__ == "__main__":
    generate_configs()