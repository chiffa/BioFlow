"""
An example of a top-level script that could be as an example of analysis of the data.
"""
from bioflow.utils.log_behavior import clear_logs
from bioflow.utils.source_dbs_download import pull_online_dbs
from bioflow.configs.main_configs import sources_location
from bioflow.db_importers.import_main import destroy_db, build_db
from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as \
    knowledge_analysis, ref_param_set, _filter
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians
import os
from bioflow.utils.smtp_log_behavior import get_smtp_logger, started_process, \
    successfully_completed, smtp_error_bail_out
from bioflow.utils.log_behavior import get_logger


log_smtp = get_smtp_logger('main_smtp_logger')

log = get_logger(__name__)

if __name__ == "__main__":
    try:
        started_process()
    except Exception as e:
        smtp_error_bail_out()
        raise e
    try:
        pass  # for syntactic reasons

        # first, let's clear logs:
        # clear_logs()

        # # if needed, clear the mongodb:
        # from bioflow.main_configs import annotome_rand_samp, interactome_rand_samp_db
        # annotome_rand_samp.drop()
        # interactome_rand_samp_db.drop()

        # # pulling the online databases
        # pull_online_dbs()

        ##########################################
        # After you've changed folders/sources above, you need to re-start python to force
        # main_configs update
        ##########################################

        # # clearing the database, if required
        # destroy_db()
        #
        # # building the neo4j database
        # build_db()

        background_bulbs_ids = []

        # Map the bulbs we are seeking to analyze

        # hits_ids, background_ids = map_and_save_gene_ids(
        #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Matrigel vs Collagen/Matrigel_vs_collagen-tumor.tsv',
        #     '')

        # hits_ids, background_ids = map_and_save_gene_ids(
        #       '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/Hits.csv',
        # #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/top_100_hum.csv',
        #      '')

        # hits_ids, background_ids = map_and_save_gene_ids(
        #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/both_HUM.csv',
        #     ''
        #     # '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/All_genes.csv'
        # )

        # hits_ids, background_bulbs_ids = map_and_save_gene_ids(
        #     'yeast_test_gene_set-glycogen_biosynthesis.tsv',
        #     '')

        # # get the bulbs ids for the nodes we would like to analyze
        # hits_ids = get_source_bulbs_ids()


        # background_bulbs_ids = get_background_bulbs_ids()

        # rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)

        # # perform the interactome analysis
        # interactome_analysis([hits_ids],
        #                      desired_depth=30,
        #                      processors=3,
        #                      background_list=background_bulbs_ids,
        #                      skip_sampling=False,
        #                      from_memoization=False)
        #
        # # perform the knowledge analysis
        # knowledge_analysis([hits_ids],
        #                    desired_depth=20,
        #                    processors=3,
        #                    param_set=ref_param_set,
        #                    skip_sampling=False)
        #
        #
        # raise Exception('debugging')

        chromosomes_directory = "//localhome//kucharav//Projects//BioFlow paper//yeast_chr_genes"
        background_file = os.path.join(chromosomes_directory, "all_genes.tab")

        # rebuild_the_laplacians()

        # perform the interactome analysis

        # currnetly we are having an issue where the name mapping generate a list that is buffered
        # => Not that much of a problem actually

        # interactome_analysis([hits_ids],
        #                      desired_depth=20,
        #                      processors=6,
        #                      background_list=background_bulbs_ids,
        #                      skip_sampling=False,
        #                      from_memoization=False)
        #
        # # perform the knowledge analysis
        # knowledge_analysis([hits_ids],
        #                    desired_depth=20,
        #                    processors=6,
        #                    param_set=ref_param_set,
        #                    skip_sampling=False)
        #
        #
        # raise Exception('debugging')

        hits_ids, background_bulbs_ids = map_and_save_gene_ids(
            'yeast_test_gene_set-glycogen_biosynthesis_w.tsv',
            '')

        log.info('debug: hits_ids parse: %s' % hits_ids)

        interactome_analysis(source_list=[hits_ids],
                     # output_destinations_list=['chr_%s' % filename[:-4]],
                     desired_depth=5,
                     processors=1,
                     background_list=background_bulbs_ids,
                     skip_sampling=False
                     )

        raise Exception('Debug Exception')


        for filename in os.listdir(chromosomes_directory):
            if filename != "all_genes.tab":

                target_file = os.path.join(chromosomes_directory, filename)
                hits_ids, background_bulbs_ids = map_and_save_gene_ids(target_file, background_file)

                # if not background_set:
                #     rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)
                #     background_set = True

                # # perform the interactome analysis

                interactome_analysis(source_list=[hits_ids[:20]],
                                     output_destinations_list=['chr_%s' % filename[:-4]],
                                     desired_depth=5,
                                     processors=1,
                                     background_list=background_bulbs_ids,
                                     skip_sampling=False
                                     )

                raise Exception('Debug Exception')

                # # perform the knowledge analysis
                knowledge_analysis(source_list=[hits_ids[20:]],
                                   output_destinations_list=['chr_%s' % filename[:-4]],
                                   desired_depth=5,
                                   processors=1,
                                   background_list=background_bulbs_ids,
                                   skip_sampling=False,
                                   )

                # raise Exception('debug')
    except Exception as e:
        try:
            log_smtp.exception(e)
        except Exception as e:
            smtp_error_bail_out()
            raise e
        raise e

    else:
        try:
            successfully_completed()
        except Exception as e:
            smtp_error_bail_out()
            raise e
