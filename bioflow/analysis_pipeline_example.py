"""
Top-Level scripts, examples of analysis pipelines.
"""

from bioflow.annotation_network.knowledge_access_analysis \
    import auto_analyze as knowledge_analysis
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.utils.io_routines import get_source_bulbs_ids
from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians

if __name__ == "__main__":
    # # first, let's clear logs:
    # clear_logs()

    # # if needed, clear the mongodb:
    # annotome_rand_samp.drop()
    # interactome_rand_samp_db.drop()

    # # setting static folders and urls for the databases
    # set_folders('/home/andrei/support')
    # # pulling the online databases
    # pull_online_dbs()

    # # setting the organism to XXXX
    # build_source_config('human')
    # raise Exception('planned interrupt')

    ##########################################
    # After you've changed folders/sources above, you need to re-start python to force
    # main_configs update
    ##########################################

    # # # clearing the database, if required
    # destroy_db()

    # # building the neo4j database
    # build_db()

    background_bulbs_ids = []

    # hits_ids, background_ids = map_and_save_gene_ids('/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Matrigel vs Collagen/Matrigel_vs_collagen-tumor.tsv',
    #                                                  '')

    # hits_ids, background_ids = map_and_save_gene_ids(
    #       '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/Hits.csv',
    # #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/top_100_hum.csv',
    #      '')

    hits_ids, background_ids = map_and_save_gene_ids(
        '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/both_HUM.csv',
        ''
        # '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/All_genes.csv'
    )

    # get the bulbs ids if the nodes we would like to analyze
    hits_ids = get_source_bulbs_ids()

    # background_bulbs_ids = get_background_bulbs_ids()

    rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)

    # perform the interactome analysis
    interactome_analysis([hits_ids],
                         desired_depth=9,
                         processors=3,
                         background_list=background_bulbs_ids,
                         skip_sampling=False,
                         from_memoization=False)

    # perform the knowledge analysis
    knowledge_analysis([hits_ids],
                       desired_depth=20,
                       processors=3,
                       param_set=ref_param_set,
                       skip_sampling=False)
