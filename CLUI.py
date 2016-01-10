__author__ = 'ank'

import click

from bioflow.configs_manager import StructureGenerator, set_folders
from bioflow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface as AnnotomeInterface, get_background
from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as knowledge_analysis
from bioflow.db_importers.import_main import build_db, destroy_db
from bioflow.main_configs import neo4j_server, annotome_rand_samp, interactome_rand_samp
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface as InteractomeInterface
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set


@click.group()
def main():
    pass


@click.command()
@click.option('--path', prompt='Please provide the path where the data will be stored')
@click.option('--neo4jserver', default='http://localhost:7474',
              prompt='Please provide the address and port on which neo4j is currently serving')
@click.option('--mongoserver', default='mongodb://localhost:27017/',
              prompt='Please provide the address and port on which mongodb is currently serving')
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
@click.confirmation_option(help='Are you sure you want to purge this neo4j database instance?'
                                ' You will have to re-import all the data for this to work properly')
def purgeneo4j():
    """
    Wipes the neo4j organism-specific database

    :return:
    """
    print 'neo4j will start purging the master database. It will take some time to finish.' \
          ' Please do not close the shell. You can supervise the progress through' \
          ' %s/webadmin interface' % neo4j_server
    destroy_db()

@click.command()
@click.confirmation_option(help='Are you sure you want to start loading the neo4j '
                                'database? The process might take several hours or days')
def loadneo4j():
    """
    Loads the information from external database into the master repositry inside neo4j

    :return:
    """
    print 'neo4j will start loading data into the master database. It will take a couple ' \
          'of hours to finish. Please do not close the shell. You can supervise the progress' \
          ' through %s/webadmin interface' % neo4j_server
    build_db()


@click.command()
@click.option('--interactome', 'matrixtype', flag_value='interactome',
              default=True)
@click.option('--annotmap', 'matrixtype', flag_value='annotome')
def extractmatrix(matrixtype):
    if matrixtype == 'interactome':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.full_rebuild()
        print local_matrix.laplacian_matrix, local_matrix.matrix_index_2_bulbs_id
    if matrixtype == 'annotome':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.full_rebuild()
        filtr = ['biological_process']
        annot_matrix = AnnotomeInterface(filtr, local_matrix.all_uniprots_bulbs_id_list, (1, 1), True, 3)
        annot_matrix.rebuild()
        print annot_matrix.laplacian_matrix, annot_matrix.Num2GO


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
@click.argument('bioflow') ##defaults
def analyze(matrixtype, background, source, depth, processors,):
    if matrixtype == 'interactome':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.full_rebuild()
        if background is not None:
            background = get_background(background)
        # TODO: rename this function and use only one of it everywhere
        source_set = get_background(source)
        interactome_analysis(source_set, depth, processors, background)
    elif matrixtype == 'annotmap':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.full_rebuild()
        filtr = ['biological_process']
        if background is None:
            annot_matrix = AnnotomeInterface(filtr, local_matrix.all_uniprots_bulbs_id_list, (1, 1), True, 3)
        else:
            background_set = get_background(background)
            annot_matrix = AnnotomeInterface(filtr, background_set, (1, 1), True, 3)
        annot_matrix.rebuild()
        annot_matrix.store()
        # TODO: rename this function and use only one of it everywhere
        source_set = get_background(source)
        knowledge_analysis(source=source_set, KG_object=annot_matrix,
                           desired_depth=depth, processors=processors)
    print "analsysis is finished, current results are stored " \
          "in the $PROJECT_HOME/bioflow/outputs directory"


@click.command()
@click.option('--drop_type', type=click.Choice(['all', 'interactome', 'annotome']))
def purgemongo(drop_type):
    if drop_type == 'all':
        annotome_rand_samp.drop()
        interactome_rand_samp.drop()
    elif drop_type == 'interactome':
        interactome_rand_samp.drop()
    elif drop_type == 'annotome':
        annotome_rand_samp.drop()


main.add_command(initialize)
main.add_command(setorgconfs)
main.add_command(downloaddbs)
main.add_command(purgeneo4j)
main.add_command(loadneo4j)
main.add_command(extractmatrix)
main.add_command(mapids)
main.add_command(analyze)
main.add_command(purgemongo)

if __name__ == '__main__':
    main()
