"""
Top-Level scripts, examples of analysis pipelines.
"""
import os
from os import path

from bioflow.configs_manager import set_folders, build_source_config, \
    pull_online_dbs
from bioflow.annotation_network.BioKnowledgeInterface \
    import GeneOntologyInterface as AnnotomeInterface
from bioflow.annotation_network.knowledge_access_analysis \
    import auto_analyze as knowledge_analysis
from bioflow.db_importers.import_main import build_db, destroy_db
from bioflow.main_configs import neo4j_server, annotome_rand_samp, interactome_rand_samp, \
    analysis_protein_ids_csv
from bioflow.molecular_network.InteractomeInterface \
    import InteractomeInterface as InteractomeInterface
from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set, \
    cast_analysis_set_to_bulbs_ids, cast_background_set_to_bulbs_id
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.log_behavior import clear_logs


def map_and_save_gene_ids(hit_genes_location, all_detectable_genes_location=''):
    cast_analysis_set_to_bulbs_ids(hit_genes_location)
    hit_genes_ids = get_source_bulbs_ids()

    if all_detectable_genes_location:
        cast_background_set_to_bulbs_id(
            background_set_csv_location=all_detectable_genes_location,
            analysis_set_csv_location=hit_genes_location)

        all_detectable_genes_ids = get_background_bulbs_ids()

    else:
        all_detectable_genes_ids = []

    return hit_genes_ids, all_detectable_genes_ids


def rebuild_the_laplacians(all_detectable_genes=[]):

    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=False)
    local_matrix.full_rebuild()

    _filter = ['biological_process']
    ref_param_set = [_filter, all_detectable_genes, (1, 1), True, 3]  # TODO: all_detectable_genes should not be here either

    annot_matrix = AnnotomeInterface(*ref_param_set)
    annot_matrix.full_rebuild()


if __name__ == "__main__":
    # # first, let's clear logs:
    # clear_logs()

    # # if needed, clear the mongodb:
    # annotome_rand_samp.drop()
    # interactome_rand_samp.drop()

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

    hits_ids, background_ids = map_and_save_gene_ids('/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Matrigel vs Collagen/Matrigel_vs_collagen-tumor.tsv',
                                                     '')

    # hits_ids, background_ids = map_and_save_gene_ids(
    #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/Hits.csv',
    #     '')

    # hits_ids, background_ids = map_and_save_gene_ids(
    #     '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/both_HUM.csv',
    #     ''
    #     # '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TWIST1_ECAD/All_genes.csv'
    #     )

    # get the bulbs ids if the nodes we would like to analyze
    hits_ids = get_source_bulbs_ids()

    # background_bulbs_ids = get_background_bulbs_ids()

    _filter = ['biological_process']

    # rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)

    # perform the interactome analysis
    interactome_analysis([hits_ids],
                         desired_depth=12,
                         processors=4,
                         background_list=background_bulbs_ids,
                         skip_sampling=False,
                         from_memoization=False)

    # # perform the knowledge analysis
    # knowledge_analysis([hits_ids],
    #                    desired_depth=20,
    #                    processors=3,
    #                    param_set=ref_param_set,
    #                    skip_sampling=False)
