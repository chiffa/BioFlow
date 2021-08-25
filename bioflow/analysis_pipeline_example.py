"""
An example of a top-level script that could be as an example of analysis of the data.
"""
from bioflow.utils.log_behavior import clear_logs
from bioflow.utils.source_dbs_download import pull_online_dbs
from bioflow.configs.main_configs import sources_location, smtp_logging
from bioflow.db_importers.import_main import destroy_db, build_db
from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as \
    knowledge_analysis
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians, \
    generate_random_weights
import os
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


if smtp_logging:
    from bioflow.utils.smtp_log_behavior import mail_handler
    log.addHandler(mail_handler)


if __name__ == "__main__":

    # first, let's clear logs:
    # clear_logs()

    # # if needed, clear the mongodb:
    # from bioflow.main_configs import annotome_rand_samp, interactome_rand_samp_db
    # annotome_rand_samp.drop()
    # interactome_rand_samp_db.drop()

    # # pulling the online databases
    # pull_online_dbs()


    # # clearing the database, if required
    # destroy_db()
    #
    # # building the neo4j database
    # build_db()

    background_internal_ids = []

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

    # hits_ids, background_internal_ids = map_and_save_gene_ids(
    #     'yeast_test_gene_set-glycogen_biosynthesis.tsv',
    #     '')

    # # get the bulbs ids for the nodes we would like to analyze
    # hits_ids = get_source_bulbs_ids()


    # background_internal_ids = get_background_bulbs_ids()

    # rebuild_the_laplacians(all_detectable_genes=background_internal_ids)

    # # perform the interactome analysis
    # interactome_analysis([hits_ids],
    #                      random_samples_to_test_against=30,
    #                      processors=3,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False,
    #                      from_memoization=False)
    #
    # # perform the knowledge analysis
    # knowledge_analysis([hits_ids],
    #                    random_samples_to_test_against=20,
    #                    processors=3,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
    #
    #
    # raise Exception('debugging')

    # generate_random_weights(os.path.join(chromosomes_directory, "all_genes.tab"),
    #                         'yeast_test_weighted_background.tsv')

    # rebuild_the_laplacians()

    # perform the interactome analysis

    # currnetly we are having an issue where the name mapping generate a list that is buffered
    # => Not that much of a problem actually

    # interactome_analysis([hits_ids],
    #                      random_samples_to_test_against=20,
    #                      processors=6,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False,
    #                      from_memoization=False)
    #
    # # perform the knowledge analysis
    # knowledge_analysis([hits_ids],
    #                    random_samples_to_test_against=20,
    #                    processors=6,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
    #
    #
    # raise Exception('debugging')

    # hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
    #     'yeast_test_gene_set-glycogen_biosynthesis_w.tsv',
    #     '')
    #
    # log.debug('debug: hits_ids parse: %s' % hits_ids)
    #
    # interactome_analysis(source_list=hits_ids,
    #                      secondary_source_list=sec_hit_ids,
    #                      # output_destinations_list=['chr_%s' % filename[:-4]],
    #                      random_samples_to_test_against=5,
    #                      processors=1,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False
    #                      )
    #
    # raise Exception('Debug Exception')

    # hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
    #     ('yeast_test_gene_set-glycogen_biosynthesis_tsw_1.tsv',
    #      'yeast_test_gene_set-glycogen_biosynthesis_tsw_2.tsv'),
    #      'yeast_test_weighted_background.tsv')

    test_folder_root = '/localhome/kucharav/Ewald Lab validation/distilled/'

    matrigel_collagen_tumor = 'Collagen-Matrigel/Matrigel_vs_collagen-tumor.tsv'
    matrigel_collagen_normal = 'Collagen-Matrigel/Matrigel_vs_collagen-normal.tsv'

    kp_km_mouse = 'Kp-Km/top_100.csv'
    kp_km_human = 'Kp-Km/top_100_hum.csv'
    kp_km_20_rd_drop = 'Kp-Km/top_80 - 20p drop-out.csv'

    K14_base = 'K14/All_genes.csv'
    K14_human = 'K14/both_HUM.csv'
    K14_mouse = 'K14/both.csv'

    TWIST_1_base = 'TWIST-1/All_genes.csv'
    TWIST_1_hits = 'TWIST-1/Hits.csv'


    second_folder_root = '/localhome/kucharav/bioflow validation data'

    breast_cancer_1q = 'breast cancer aneuploidy/1q.txt'
    breast_cancer_8p = 'breast cancer aneuploidy/8p.txt'
    breast_cancer_8q = 'breast cancer aneuploidy/8q.txt'
    breast_cancer_8q2 = 'breast cancer aneuploidy/8q - second half.txt'
    breast_cancer_18q = 'breast cancer aneuploidy/18q.txt'
    breast_cancer_20q = 'breast cancer aneuploidy/20q.txt'

    Linhao_no_import_hits = '"yeast -  Linhao no import screens"/hits.csv'
    Linhao_no_import_background = '"yeast -  Linhao no import screens"background.csv'

    Linhao_MudPIT_hits = 'yeast - Linhao MudPIT/background.csv'
    Linhao_MudPIT_background = 'yeast - Linhao MudPIT/hits_WT-MUT_vs_CTRL.csv'


    active_load = os.path.join(second_folder_root, breast_cancer_20q)
    active_back = ''
    # active_back = os.path.join(second_folder_root, TWIST_1_base)

    hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
         active_load, active_back)

    # log.debug('debug: hits_ids parse: %s' % hits_ids)

    log.info('Registration: run from %s' % active_load)

    interactome_analysis(source_list=hits_ids,
                         secondary_source_list=sec_hit_ids,
                         output_destinations_list=['BR_C_8q'],
                         random_samples_to_test_against=25,
                         processors=1,
                         background_list=background_internal_ids,
                         skip_sampling=False
                         )

    knowledge_analysis(source_list=hits_ids,
                       secondary_source_list=sec_hit_ids,
                       output_destinations_list=['BR_C_8q'],
                       random_samples_to_test_against=25,
                       processors=1,
                       background_list=background_internal_ids,
                       skip_sampling=False,
                       )

    # chromosomes_directory = "//localhome//kucharav//Projects//BioFlow paper//yeast_chr_genes"
    # background_file = os.path.join(chromosomes_directory, "all_genes.tab")
    #
    # chromosomes_files = []
    # name_files = []
    # for filename in os.listdir(chromosomes_directory):
    #     if filename != "all_genes.tab":
    #         chromosomes_files.append(os.path.join(chromosomes_directory, filename))
    #         name_files.append('chr_%s' % filename[:-4])
    #
    # hits_ids, sec_hit_ids, background_ids = map_and_save_gene_ids(chromosomes_files, background_file)
    #
    # interactome_analysis(source_list=hits_ids,
    #                      output_destinations_list=name_files,
    #                      random_samples_to_test_against=25,
    #                      processors=1,
    #                      background_list=background_ids,
    #                      skip_sampling=False
    #                      )
    #
    # knowledge_analysis(source_list=hits_ids,
    #                    output_destinations_list=name_files,
    #                    random_samples_to_test_against=25,
    #                    processors=1,
    #                    background_list=background_ids,
    #                    skip_sampling=False,
    #                    )


    # # compile chromosome file names
    #
    # for filename in os.listdir(chromosomes_directory):
    #     if filename != "all_genes.tab":
    #
    #         target_file = os.path.join(chromosomes_directory, filename)
    #         hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(target_file,
    #                                                                                background_file)
    #
    #         # if not background_set:
    #         #     rebuild_the_laplacians(all_detectable_genes=background_internal_ids)
    #         #     background_set = True
    #
    #         # # perform the interactome analysis
    #
    #         interactome_analysis(source_list=hits_ids,
    #                              output_destinations_list=['chr_%s' % filename[:-4]],
    #                              random_samples_to_test_against=5,
    #                              processors=1,
    #                              background_list=background_internal_ids,
    #                              skip_sampling=False
    #                              )
    #
    #         # # perform the knowledge analysis
    #         knowledge_analysis(source_list=hits_ids,
    #                            output_destinations_list=['chr_%s' % filename[:-4]],
    #                            random_samples_to_test_against=5,
    #                            processors=1,
    #                            background_list=background_internal_ids,
    #                            skip_sampling=False,
    #                            )
    #
    #         # raise Exception('debug')
