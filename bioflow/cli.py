"""
This modules manages the command line interface
"""
import click


# CURRENTPASS: add smtp logging if enabled

def print_version(ctx, value):
    from bioflow import __version__
    if not value or ctx.resilient_parsing:
        return
    click.echo(__version__)
    ctx.exit()


@click.command()
def about():
    """
    Information about the project

    :return:
    """
    from bioflow import __version__, __author__, __author_mail__, __current_year__

    click.echo('BioFlow \n'
               '\tversion: %s \n'
               '\tauthor: %s \n'
               '\tcontact: %s \n'
               '\tLicense: %s\n'
               '\tcite: %s\n'
               '\n 2013-%s Andrei Kucharavy all rights reserved' %
               (__version__, __author__, __author_mail__, 'BSD 3-clause',
                'https://github.com/chiffa/BioFlow', __current_year__))


@click.group()
@click.option('--version', '-v', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def main():
    pass


@click.command()
@click.confirmation_option(help='Pulling online databases. '
                                'Please make sure you initialized the project and you are ready'
                                'to wait for a while for hte download to complete. '
                                'Some files are large (up to 3 Gb).'
                                'You can perform this step manually (cf documentation).')
def downloaddbs():
    """
    Downloads the databases automatically
    \f

    :return:
    """
    from bioflow.utils.source_dbs_download import pull_online_dbs
    pull_online_dbs()


@click.command()
@click.confirmation_option(help='Are you sure you want to purge this neo4j database instance?'
                                ' You will have to re-import all the data'
                                'for this to work properly')
def purgeneo4j():
    """
    Wipes the neo4j organism-specific database
    \f

    :return:
    """
    print('neo4j will start purging the master database. It will take some time to finish.' \
          ' Please do not close the shell')
    from bioflow.db_importers.import_main import destroy_db
    destroy_db()


@click.command()
@click.confirmation_option(help='Are you sure you want to start loading the neo4j '
                                'database? The process might take several hours or days')
def loadneo4j():
    """
    Loads the information from external database into the main knowledge repository inside neo4j
    \f

    :return:
    """
    print('neo4j will start loading data into the master database. It will take a couple ' \
          'of hours to finish. Please do not close the shell.')
    from bioflow.db_importers.import_main import build_db
    build_db()


@click.command()
@click.option('--background', default='', help='path to file of IDs of all genes detectable by a method')
@click.argument('source')
def mapsource(source, background):
    """
    Sets the source and background files that will be uses in the analysis.

    The argument source is a path to a file containing the IDs of all genes considered as a hit

    Preferred formats
    are HGCN gene names (TP53), Uniprot gene names (P53_HUMAN) or Uniprot Accession numbers (
    P04637).
    Other sources, such as ENSEMBL or PDB IDs are supported as well
    \f

    :param source:
    :param background:
    :return:
    """
    from bioflow.utils.top_level import map_and_save_gene_ids
    map_and_save_gene_ids(source, background)


@click.command()
def rebuildlaplacians():
    """
    Extracts the Laplacian matrices from the master graph database.
    \f

    :return:
    """
    from bioflow.utils.top_level import rebuild_the_laplacians
    rebuild_the_laplacians()


@click.command()
@click.option('--collection', type=click.Choice(['all', 'interactome', 'annotome']), default='all')
def purgemongo(collection):
    """
    purges the mongodb collection currently used to store all the information.
    \f

    :param collection:
    :return:
    """
    from bioflow.sample_storage.mongodb import drop_all_interactome_rand_samp
    from bioflow.sample_storage.mongodb import drop_all_annotome_rand_samp

    if collection == 'all':
        drop_all_annotome_rand_samp()
        drop_all_interactome_rand_samp()
    elif collection == 'interactome':
        drop_all_interactome_rand_samp()
    elif collection == 'annotome':
        drop_all_annotome_rand_samp()


@click.command()
@click.option('--matrix', type=click.Choice(['all', 'interactome', 'annotome']), default='all',
              help='analyse molecular entities alone (interactome), annotation entities alone ('
                   'annotome) or both')
@click.option('--depth', default=25, help='random samples used to infer flow pattern significance')
@click.option('--processors', default=1, help='processor cores used in flow patterns calculation')
@click.option('--skipsampling', default=False, help='if True, skips random sparse_sampling step')
@click.option('--background', default=False, help='if True, uses the background for sparse_sampling')
@click.option('--name', default='', help='name of the experiment')
def analyze(name, matrix, depth, processors, skipsampling, background):
    """
    Performs the analysis of the information flow

    :param name:
    :param matrix:
    :param depth:
    :param processors:
    :param skipsampling:
    :param background:
    :return:
    """
    from bioflow.utils.io_routines import get_background_bulbs_ids, get_source_bulbs_ids
    from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
    from bioflow.annotation_network.knowledge_access_analysis \
        import auto_analyze as knowledge_analysis

    source = get_source_bulbs_ids()

    if background:
        background = get_background_bulbs_ids()

    else:
        background = []

    if name:
        name = [name]
    else:
        name = []

    if matrix != 'annotome':
        # # perform the interactome analysis
        interactome_analysis([source],
                             name,
                             desired_depth=depth,
                             processors=processors,
                             background_list=background,
                             skip_sampling=skipsampling
                             )

    if matrix != 'interactome':
        # # perform the knowledge analysis
        knowledge_analysis([source],
                           name,
                           desired_depth=depth,
                           processors=processors,
                           background_list=background,
                           skip_sampling=skipsampling,
                            )




# @click.command()
# @click.option('--depth', default=24, help='random samples used to infer flow pattern significance')
# @click.option('--processors', default=3, help='processor cores used in flow patterns calculation')
# @click.option('--skipsampling', default=False, help='if True, skips random sparse_sampling step')
# @click.option('--background', default=False, help='if True, skips hits sample flow computation ')
# def interactomeanalysis(depth, processors, skipsampling, background):
#     """
#     Performs interactome analysis given background set given earlier.
#     \f
#
#     :param depth:
#     :param processors:
#     :param skipsampling:
#     :param skiphitflow:
#     :return:
#     """
#     from bioflow.utils.io_routines import get_background_bulbs_ids, get_source_bulbs_ids
#     from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
#
#     source_bulbs_ids = get_source_bulbs_ids()
#     background_bulbs_ids = get_background_bulbs_ids()
#
#     interactome_analysis([source_bulbs_ids],
#                          desired_depth=depth,
#                          processors=processors,
#                          background_list=background_bulbs_ids,
#                          skip_sampling=skipsampling,
#                          from_memoization=skiphitflow)
#
#
# @click.command()
# @click.option('--depth', default=24, help='random samples used to infer flow pattern significance')
# @click.option('--processors', default=3, help='processor cores used in flow patterns calculation')
# @click.option('--skipsampling', default=False, help='if True, skips random sparse_sampling step')
# def knowledgeanalysis(depth, processors, skipsampling):
#     """
#     Performs annotome analysis given background set given earlier.
#     \f
#
#     :param depth:
#     :param processors:
#     :param skipsampling:
#     :return:
#     """
#     from bioflow.utils.io_routines import get_source_bulbs_ids
#     from bioflow.annotation_network.knowledge_access_analysis \
#         import auto_analyze as knowledge_analysis
#
#     source_bulbs_ids = get_source_bulbs_ids()
#
#     knowledge_analysis([source_bulbs_ids],
#                        desired_depth=depth,
#                        processors=processors,
#                        skip_sampling=skipsampling)


main.add_command(downloaddbs)
main.add_command(purgeneo4j)
main.add_command(loadneo4j)
main.add_command(rebuildlaplacians)
main.add_command(purgemongo)
main.add_command(mapsource)
main.add_command(analyze)
main.add_command(about)

if __name__ == '__main__':
    main()
