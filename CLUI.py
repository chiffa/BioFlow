__author__ = 'ank'

import click

from BioFlow.configs_manager import StructureGenerator, set_folders
from BioFlow.annotation_network.BioKnowledgeInterface import GO_Interface as AnnotomeInterface, get_background
from BioFlow.annotation_network.knowledge_access_analysis import auto_analyze as knowledge_analysis
from BioFlow.db_importers.import_main import build_db, destroy_db
from BioFlow.main_configs import neo4j_server
from BioFlow.molecular_network.InteractomeInterface import MatrixGetter as InteractomeInterface
from BioFlow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from BioFlow.neo4j_db.db_io_routines import look_up_annotation_set


@click.group()
def truegird():
    pass


@click.command()
@click.option('--path', prompt='Please provide the path where the data will be stored')
@click.option('--neo4jserver', default='http://localhost:7474', prompt='Please provide the address and port on which neo4j is currently serving')
@click.option('--mongoserver', default='mongodb://localhost:27017/', prompt='Please provide the address and port on which mongodb is currently serving')
def initialize(path, neo4jserver, mongoserver):
    """
    Initialized the working environement

    :param path: path where the external data sources are expected to be stored
    :param neo4jserver: neo4j server adress and port
    :param mongoserver: mongodb server adress and port
    :return:
    """
    set_folders(path, neo4jserver, mongoserver)


@click.command()
@click.confirmation_option(help='Pulling online databases. Please make sure you initialized the project and you are ready\n'+
                                'to wait for a while for hte download to complete. Some files are large (up to 3 Gb).\n'+
                                'You can perform this step manually (cf documentation). Are you sure you want to continue?')
def downloaddbs():
    """
    Downloads the databases automatically

    :return:
    """
    StructureGenerator.pull_online_dbs()


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
@click.confirmation_option(help='Are you sure you want to purge this neo4j database instance? You will have to re-import all the data for this to work properly')
def purgeneo4j():
    """
    Wipes the neo4j organism-specific database

    :return:
    """
    print 'neo4j will start purging the master database. It will take some time to finish. Please do not close the shell. You can supervise the progress through %s/webadmin interface' % neo4j_server
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
        local_matrix = InteractomeInterface(Connexity_Aware=True, full_impact=True)
        local_matrix.full_rebuild()
        print local_matrix.Conductance_Matrix, local_matrix.MatrixNumber2NodeID
    if matrixtype == 'annotome':
        local_matrix = InteractomeInterface(Connexity_Aware=True, full_impact=True)
        local_matrix.full_rebuild()
        filtr = ['biological_process']
        annot_matrix = AnnotomeInterface(filtr, local_matrix.Uniprot_complete, (1, 1), True, 3)
        annot_matrix.rebuild()
        print annot_matrix.Laplacian_matrix, annot_matrix.Num2GO


@click.command()
@click.argument('idlist')
def mapids(idlist):
    print look_up_annotation_set(idlist)


@click.command()
@click.option('--interactome', 'matrixtype', flag_value='interactome',
              default=True)
@click.option('--annotmap', 'matrixtype', flag_value='annotmap')
@click.option('--background', default=None) ##defaults
@click.option('--depth', default=100) ##defaults
@click.option('--processors', default=2) ##defaults
@click.argument('source') ##defaults
def analyze(matrixtype, background, source, depth, processors,):
    if matrixtype == 'interactome':
        local_matrix = InteractomeInterface(Connexity_Aware=True, full_impact=True)
        local_matrix.full_rebuild()
        if background is not None:
            background = get_background(background)
        source_set = get_background(source)  # TODO: rename this function and use only one of it everywhere
        interactome_analysis(source_set, depth, processors, background)
    elif matrixtype == 'annotmap':
        local_matrix = InteractomeInterface(Connexity_Aware=True, full_impact=True)
        local_matrix.full_rebuild()
        filtr = ['biological_process']
        if background is None:
            annot_matrix = AnnotomeInterface(filtr, local_matrix.Uniprot_complete, (1, 1), True, 3)
        else:
            background_set = get_background(background)
            annot_matrix = AnnotomeInterface(filtr, background_set, (1, 1), True, 3)
        annot_matrix.rebuild()
        annot_matrix.store()
        source_set = get_background(source)  # TODO: rename this function and use only one of it everywhere
        knowledge_analysis(source=source_set, KG_object=annot_matrix, desired_depth=depth, processors=processors)
    print "analsysis is finished, current results are stored in the $PROJECT_HOME/BioFlow/outputs directory"


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