__author__ = 'ank'

import click
from PolyPharma.Utils.ConfigsIO import StructureGenerator, edit_confile
from os.path import abspath

@click.group()
def truegird():
    pass


@click.command()
@click.option('--path', prompt='Please provide the path where the data will be stored')
@click.option('--neo4jserver', default='http://localhost:7474', prompt='Please provide the adresss and port on which neo4j is currently serving')
@click.option('--mongoserver', default='mongodb://localhost:27017/', prompt='Please provide the adresss and port on which mongodb is currently serving')
def initialize(path, neo4jserver, mongoserver):
    """
    Initialized the working environement

    :param path: path where the external data sources are expected to be stored
    :param neo4jserver: neo4j server adress and port
    :param mongoserver: mongodb server adress and port
    :return:
    """
    edit_confile('servers', 'PRODUCTION', 'base_folder', abspath(path))
    edit_confile('servers', 'PRODUCTION', 'server_neo4j', neo4jserver)
    edit_confile('servers', 'PRODUCTION', 'mongodb_server', mongoserver)


@click.command()
@click.option('--organism', type=click.Choice(['mouse', 'human', 'yeast']))
def setorgconfigs(organism):
    """
    Sets organism-specific configurations

    :param organism:
    :return:
    """
    StructureGenerator.build_source_config(organism)


@click.command()
@click.confirmation_option(help='Are you sure you want to purge this neo4j database insance? You will have to re-import all the data for this to work properly')
def purgeneo4j():
    """
    Wipes the neo4j organism-specific database

    :return:
    """
    pass

@click.command()
@click.confirmation_option(help='Are you sure you want to start loading the neo4j database? The process might take several hours or days')
def loadneo4j():
    """
    Loads the information from external database into the master repositry inside neo4j

    :return:
    """
    print 'neo4j will start loading. Please do not close the shell. You can supervise progress through %s/webadmin interface' % neo4j_server


#TODO: purge mongodb



truegird.add_command(initialize)
truegird.add_command(setorgconfigs)
truegird.add_command(purgeneo4j)
truegird.add_command(loadneo4j)

if __name__ == '__main__':
    truegird()