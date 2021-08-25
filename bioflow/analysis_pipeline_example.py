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
from pathlib import Path
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


if smtp_logging:
    from bioflow.utils.smtp_log_behavior import mail_handler
    log.addHandler(mail_handler)


if __name__ == "__main__":

    # ==========================================
    # Setup
    # ==========================================

    # # first, let's clear logs:
    # clear_logs()

    # # pulling the online databases
    # pull_online_dbs()

    # # clear the main database, if needed
    # destroy_db()

    # # rebuild the main database, if needed
    # build_db()

    # # finally rebuild the laplacians:
    # rebuild_the_laplacians()

    # # in case we are constantly switch around in code between using or not background
    # background_internal_ids = []

    # ==============================================================
    # Perform integration tests on yeast glycogen biosynthesis tests
    # ===============================================================

    # ------------------------
    # just the unweighted set
    # ------------------------

    # hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
    #     'yeast_test_gene_set-glycogen_biosynthesis.tsv',
    #     '')

    # # perform the interactome analysis
    # interactome_analysis(hits_ids,
    #                      random_samples_to_test_against=20,
    #                      processors=3,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False,
    #                      from_memoization=False)
    #
    # # perform the knowledge analysis
    # knowledge_analysis(hits_ids,
    #                    random_samples_to_test_against=20,
    #                    processors=3,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
    #
    #

    # ------------------------
    # weighted set
    # ------------------------

    # hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
    #     'yeast_test_gene_set-glycogen_biosynthesis_w.tsv',
    #     '')

    # # perform the interactome analysis
    # interactome_analysis(hits_ids,
    #                      random_samples_to_test_against=20,
    #                      processors=3,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False,
    #                      from_memoization=False)
    #
    # # perform the knowledge analysis
    # knowledge_analysis(hits_ids,
    #                    random_samples_to_test_against=20,
    #                    processors=3,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
    #
    #

    # ---------------------------------
    # secondary and background weighted
    # ---------------------------------

    # hits_ids, sec_hit_ids, background_internal_ids = map_and_save_gene_ids(
    #     ('yeast_test_gene_set-glycogen_biosynthesis_tsw_1.tsv',
    #      'yeast_test_gene_set-glycogen_biosynthesis_tsw_2.tsv'),
    #      'yeast_test_weighted_background.tsv')

    # # perform the interactome analysis
    # interactome_analysis(hits_ids,
    #                      random_samples_to_test_against=20,
    #                      processors=3,
    #                      background_list=background_internal_ids,
    #                      skip_sampling=False,
    #                      from_memoization=False)
    #
    # # perform the knowledge analysis
    # knowledge_analysis(hits_ids,
    #                    random_samples_to_test_against=20,
    #                    processors=3,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
    #
    #

    # ===============================================================
    # Legacy Experiments
    # ===============================================================

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

    # ===============================================================
    # Ablation Experiments
    # ===============================================================
    root_directory = "C:\\Users\\Andrei\Dropbox\\workspaces\\JHU\\" \
                            "Ewald Lab\\Kp_Km data\\"

    root_dir_path = Path(root_directory)

    ablation_experiments_rail = ['lowest_5_percent_removed', 'lowest_5_percent_set_to_random',
                                 'random_5_percent_removed', 'random_5_percent_set_to_random',
                                 'lowest_10_percent_removed', 'lowest_10_percent_set_to_random',
                                 'random_10_percent_removed', 'random_10_percent_set_to_random',
                                 'lowest_20_percent_removed', 'lowest_20_percent_set_to_random',
                                 'random_20_percent_removed', 'random_20_percent_set_to_random',
                                 'lowest_50_percent_removed', 'lowest_50_percent_set_to_random',
                                 'random_50_percent_removed', 'random_50_percent_set_to_random',
                                 'no_weights',
                                 'no_weights_lowest_5_percent_removed',
                                 'no_weights_lowest_5_percent_set_to_random',
                                 'no_weights_random_5_percent_removed',
                                 'no_weights_random_5_percent_set_to_random',
                                 'no_weights_lowest_10_percent_removed',
                                 'no_weights_lowest_10_percent_set_to_random',
                                 'no_weights_random_10_percent_removed',
                                 'no_weights_random_10_percent_set_to_random',
                                 'no_weights_lowest_20_percent_removed'
                                 'no_weights_lowest_20_percent_set_to_random',
                                 'no_weights_random_20_percent_removed'
                                 'no_weights_random_20_percent_set_to_random',
                                 'no_weights_lowest_50_percent_removed',
                                 'no_weights_lowest_50_percent_set_to_random',
                                 'no_weights_random_50_percent_removed',
                                 'no_weights_random_50_percent_set_to_random']

    background_file = root_dir_path.join('mouse_genes_background.txt')
    reference_file = root_dir_path.join('mouse_weighted_abs_log-fold.txt')

    fname = reference_file.stem
    fname = fname + '_ablations'
    storage_folder = root_dir_path.joinpath(fname)

    path_files = [storage_folder.joinpath(stem + '.tsv') for stem in ablation_experiments_rail]

    hits_ids, sec_hit_ids, background_ids = map_and_save_gene_ids(path_files,
                                                                  background_file)

    interactome_analysis(source_list=hits_ids,
                         output_destinations_list=ablation_experiments_rail,
                         random_samples_to_test_against=25,
                         processors=1,
                         background_list=background_ids,
                         skip_sampling=False
                         )

    knowledge_analysis(source_list=hits_ids,
                       output_destinations_list=ablation_experiments_rail,
                       random_samples_to_test_against=25,
                       processors=1,
                       background_list=background_ids,
                       skip_sampling=False,
                       )


    #===============================================================
    # Chromosome Experiments
    #===============================================================

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
