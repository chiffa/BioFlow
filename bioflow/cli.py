"""
Contains the access to the command line interface of the application.
"""
import click
from bioflow.configs_manager import set_folders, build_source_config, \
    pull_online_dbs
from bioflow.annotation_network.BioKnowledgeInterface \
    import GeneOntologyInterface as AnnotomeInterface
from bioflow.utils.io_routines import get_background_bulbs_ids, get_source_bulbs_ids
from bioflow.annotation_network.knowledge_access_analysis \
    import auto_analyze as knowledge_analysis
from bioflow.db_importers.import_main import build_db, destroy_db
from bioflow.main_configs import neo4j_server, annotome_rand_samp, interactome_rand_samp_db
from bioflow.molecular_network.InteractomeInterface \
    import InteractomeInterface as InteractomeInterface
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set, \
    cast_analysis_set_to_bulbs_ids, cast_background_set_to_bulbs_id


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('0.2.3')
    ctx.exit()

@click.group()
@click.option('--version', '-v', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
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
@click.confirmation_option(help='Pulling online databases. '
                                'Please make sure you initialized the project and you are ready'
                                'to wait for a while for hte download to complete. '
                                'Some files are large (up to 3 Gb).'
                                'You can perform this step manually (cf documentation).'
                                'Are you sure you want to continue?')
def downloaddbs():
    """
    Downloads the databases automatically

    :return:
    """
    pull_online_dbs()


@click.command()
@click.option('--organism', type=click.Choice(['mouse', 'human', 'yeast']))
def setorg(organism):
    """
    Sets organism-specific configurations

    :param organism:
    :return:
    """
    build_source_config(organism)


@click.command()
@click.confirmation_option(help='Are you sure you want to purge this neo4j database instance?'
                                ' You will have to re-import all the data'
                                'for this to work properly')
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
@click.argument('source')
def setsource(source):
    """
    Sets the source to analyze
    :param source:
    :return:
    """
    cast_analysis_set_to_bulbs_ids(source)


@click.command()
@click.argument('background')
@click.argument('source')
def setbackground(background, source):
    """
    Sets the background that would be later used in the analysis
    :param background:
    :param source:
    :return:
    """
    cast_background_set_to_bulbs_id(background, source)


@click.command()
@click.option('--interactome', 'matrixtype', flag_value='interactome',
              default=True)
@click.option('--annotome', 'matrixtype', flag_value='annotome')
def extractmatrix(matrixtype):
    """
    Extracts the matrix interface object for the computation routine.
    :param matrixtype:
    :return:
    """
    if matrixtype == 'interactome':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.full_rebuild()

    if matrixtype == 'annotome':
        local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
        local_matrix.fast_load()
        ref_param_set = [['biological_process'], get_background_bulbs_ids(), (1, 1), True, 3]
        annot_matrix = AnnotomeInterface(*ref_param_set)
        annot_matrix.full_rebuild()


@click.command()
@click.argument('idlist')
def mapids(idlist):
    """
    Maps the identifiers from the list to the bulbs node numbers in the database

    :param idlist:
    :return:
    """
    print look_up_annotation_set(idlist)


@click.command()
@click.option('--collection', type=click.Choice(['all', 'interactome', 'annotome']), default='all')
def purgemongo(collection):
    """
    purges the mongodb collection currently used to store all the information.

    :param collection:
    :return:
    """
    if collection == 'all':
        annotome_rand_samp.drop()
        interactome_rand_samp_db.drop()
    elif collection == 'interactome':
        interactome_rand_samp_db.drop()
    elif collection == 'annotome':
        annotome_rand_samp.drop()


@click.command()
@click.option('--matrixtype', type=click.Choice(['interactome', 'annotome']))
@click.option('--depth', default=24)  # defaults
@click.option('--processors', default=3)  # defaults
def analyze(matrixtype, depth, processors,):
    """
    Performs an analysis of the type given by matrixtype with the given background set and

    :param matrixtype:
    :param depth:
    :param processors:
    :return:
    """
    source_bulbs_ids = get_source_bulbs_ids()
    background_bulbs_ids = get_background_bulbs_ids()

    # TODO: CRICIAL: inject background usage when background switch is available.
    # Refer to the analysis pipeline example for an example
    interactome_interface_instance = InteractomeInterface(main_connex_only=True, full_impact=True)
    interactome_interface_instance.fast_load()
    ref_param_set = [['biological_process'], background_bulbs_ids, (1, 1), True, 3]

    if matrixtype == 'interactome':
        interactome_analysis(source_bulbs_ids, depth, processors, background_bulbs_ids)

    elif matrixtype == 'annotome':
        knowledge_analysis(source=source_bulbs_ids,
                           desired_depth=depth,
                           processors=processors,
                           param_set=ref_param_set)

    print "analsysis is finished, current results are stored " \
          "in the outputs directory"


main.add_command(initialize)
main.add_command(setorg)
main.add_command(downloaddbs)
main.add_command(purgeneo4j)
main.add_command(loadneo4j)
main.add_command(extractmatrix)
main.add_command(mapids)
main.add_command(analyze)
main.add_command(purgemongo)
main.add_command(setsource)
main.add_command(setbackground)


if __name__ == '__main__':
    main()
