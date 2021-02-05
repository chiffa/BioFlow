"""
Top-Level scripts, examples of analysis pipelines.
"""
from bioflow.configs_manager import set_folders, build_source_config, pull_online_dbs
from bioflow.db_importers.import_main import build_db, destroy_db
from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as \
    knowledge_analysis, ref_param_set, _filter, _correlation_factors
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians
from bioflow.utils.log_behavior import clear_logs
from bioflow.user_configs import sources_location
import os


if __name__ == "__main__":
    pass  # for syntactic reasons

    # first, let's clear logs:
    # clear_logs()

    # # if needed, clear the mongodb:
    # from bioflow.main_configs import annotome_rand_samp, interactome_rand_samp_db
    # annotome_rand_samp.drop()
    # interactome_rand_samp_db.drop()

    # # setting static folders and urls for the databases
    # set_folders(sources_location)

    # # pulling the online databases
    # pull_online_dbs()

    # # setting the organism to XXXX
    # build_source_config('yeast')

    # raise Exception('planned interrupt')

    ##########################################
    # After you've changed folders/sources above, you need to re-start python to force
    # main_configs update
    ##########################################

    # # clearing the database, if required
    # destroy_db()

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


    background_bulbs_ids = get_background_bulbs_ids()

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

    # rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)

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


    background_set = False
    for filename in os.listdir(chromosomes_directory):
        if filename != "all_genes.tab":
            target_file = os.path.join(chromosomes_directory, filename)
            hits_ids, background_bulbs_ids = map_and_save_gene_ids(target_file, background_file)
            paramset_with_background = tuple([_filter, background_bulbs_ids, (1, 1), True, 3])

            # if not background_set:
            #     rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)
            #     background_set = True

            # # perform the interactome analysis
            interactome_analysis([hits_ids[:20]],
                                 ['chr_%s' % filename[:-4]],
                                 desired_depth=1,
                                 processors=1,
                                 background_list=background_bulbs_ids,
                                 skip_sampling=False
                                 )

            # # perform the knowledge analysis
            knowledge_analysis([hits_ids[:20]],
                               ['chr_%s' % filename[:-4]],
                               desired_depth=1,
                               processors=1,
                               param_set=paramset_with_background,
                               skip_sampling=False,
                               )

            raise Exception('debug')

