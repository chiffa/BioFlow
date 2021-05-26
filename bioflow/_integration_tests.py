from bioflow.utils.log_behavior import clear_logs
from bioflow.utils.source_dbs_download import pull_online_dbs
from bioflow.configs.main_configs import sources_location
from bioflow.db_importers.import_main import destroy_db, build_db
from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as \
    knowledge_analysis
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians, \
    generate_random_weights
import os
from bioflow.utils.smtp_log_behavior import get_smtp_logger, started_process, \
    successfully_completed, smtp_error_bail_out
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


# TODO: transform them into actual unittests that can be ran and timed
def analysis_loop_test(primary, secondary=None, background=''):

    if secondary is None:
        hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
            primary, background)

    else:
        hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
            (primary, secondary), background)


    interactome_analysis(source_list=hits_ids,
                         secondary_source_list=sec_hit_ids,
                         # output_destinations_list=['chr_%s' % filename[:-4]],
                         random_samples_to_test_against=5,
                         processors=1,
                         background_list=background_internal_ids,
                         skip_sampling=False
                         )

    knowledge_analysis(source_list=hits_ids,
                       secondary_source_list=sec_hit_ids,
                       # output_destinations_list=['chr_%s' % filename[:-4]],
                       random_samples_to_test_against=5,
                       processors=1,
                       background_list=background_internal_ids,
                       skip_sampling=False,
                       )


if __name__ == "__main__":

    # TODO: add a set org to yeast

    # pulling the online databases
    pull_online_dbs()

    # building the neo4j database
    build_db()

    # primary alone
    analysis_loop_test('yeast_test_gene_set-glycogen_biosynthesis.tsv')
    # primary with background
    analysis_loop_test('yeast_test_gene_set-glycogen_biosynthesis.tsv',
                       'test_background.tsv')
    # primary and secondary
    analysis_loop_test(('yeast_test_gene_set-glycogen_biosynthesis_ts_1.tsv',
                        'yeast_test_gene_set-glycogen_biosynthesis_ts_2.tsv'))
    # weighted primary
    analysis_loop_test('yeast_test_gene_set-glycogen_biosynthesis_w.tsv')
    # weighted primary and secondary
    analysis_loop_test(('yeast_test_gene_set-glycogen_biosynthesis_tsw_1.tsv',
                        'yeast_test_gene_set-glycogen_biosynthesis_tsw_2.tsv'))
    # weighted primary and secondary with background
    analysis_loop_test(('yeast_test_gene_set-glycogen_biosynthesis_tsw_1.tsv',
                        'yeast_test_gene_set-glycogen_biosynthesis_tsw_2.tsv'),
                       'test_weighted_background.tsv')
