__author__ = 'ank'

import click
from PolyPharma.Utils.ConfigsIO import StructureGenerator, edit_confile
from PolyPharma.configs2 import neo4j_server
from PolyPharma.neo4j_Importers.Import_commander import build_db, destroy_db
from os.path import abspath, expanduser


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
    if path[0] == r'~':
        path = expanduser(path)
    print 'setting external DBs filestore in %s' % abspath(path)
    edit_confile('servers', 'PRODUCTION', 'base_folder', abspath(path))
    edit_confile('servers', 'PRODUCTION', 'server_neo4j', neo4jserver)
    edit_confile('servers', 'PRODUCTION', 'mongodb_server', mongoserver)
    print 're-writing the organism-specific configs for the following organism: %s' % 'yeast'
    StructureGenerator.build_source_config('yeast')


@click.command()
@click.confirmation_option(help='Pulling online databases. Please make sure you initialized the project and you are ready\n'+
                                'to wait for a while for hte download to complete. Some files are large (up to 3 Gb).\n'+
                                'You can perform this step manually (cf documetnation). Are you sure you want to continue?')
def downloaddbs():
    """
    Downloads the databases automatically

    :return:
    """
    StructureGenerator.pull_online_DBs()


@click.command()
@click.option('--organism', type=click.Choice(['mouse', 'human', 'yeast']))
def setorgconfs(organism):
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
    print 'neo4j will start purging the master databse. It will take some time to finish. Please do not close the shell. You can supervise the progress through %s/webadmin interface' % neo4j_server
    destroy_db()

@click.command()
@click.confirmation_option(help='Are you sure you want to start loading the neo4j database? The process might take several hours or days')
def loadneo4j():
    """
    Loads the information from external database into the master repositry inside neo4j

    :return:
    """
    print 'neo4j will start loading data into the master database. It will take a couple of hours to finish. Please do not close the shell. You can supervise the progress through %s/webadmin interface' % neo4j_server
    build_db()


@click.command()
@click.option('--interactome', 'matrixtype', flag_value='interactome',
              default=True)
@click.option('--annotmap', 'matrixtype', flag_value='annotome')
def extractmatrix(matrixtype):
    if matrixtype == 'interactome':
        pass
    if matrixtype == 'annotome':
        pass


@click.command()
@click.argument('idlist')
def mapids(idlist):
    pass


@click.command()
@click.option('--interactome', 'matrixtype', flag_value='interactome',
              default=True)
@click.option('--annotmap', 'matrixtype', flag_value='annotome')
@click.option('--background') ##defaults
@click.option('--depth') ##defaults
@click.option('--processors') ##defaults
@click.argument('idlist') ##defaults
def analyze(matrixtype, background, idlist, depth, processors,):
    pass


#TODO: add purge mongodb operation

truegird.add_command(initialize)
truegird.add_command(setorgconfs)
truegird.add_command(downloaddbs)
truegird.add_command(purgeneo4j)
truegird.add_command(loadneo4j)
truegird.add_command(extractmatrix)
truegird.add_command(mapids)
truegird.add_command(analyze)

if __name__ == '__main__':
    truegird()