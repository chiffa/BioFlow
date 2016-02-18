"""
Top-Level scripts, examples of analysis pipelines.
"""
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


if __name__ == "__main__":
    # # first, let's clear logs:
    # clear_logs()
    #
    # # if needed, clear the mongodb:
    # annotome_rand_samp.drop()
    # interactome_rand_samp.drop()

    # # setting static folders and urls for the databases
    # set_folders('/home/andrei/support')
    # # pulling the online databases
    # pull_online_dbs()
    # # setting the organism to yeast
    # build_source_config('human')

    # # clearing the database, if required
    # destroy_db()

    # # building the neo4j database
    # build_db()

    # # set the source file of the ids of perturbed proteins and background set:
    # cast_analysis_set_to_bulbs_ids("/home/andrei/HS_30_Linhao_outliers.txt")
    # cast_background_set_to_bulbs_id(
    #     background_set_csv_location=None,
    #     analysis_set_csv_location="/home/andrei/HS_30_Linhao_outliers.txt")

    # get the bulbs ids oif the nodes we would like to analyze
    source_bulbs_ids = get_source_bulbs_ids()
    background_bulbs_ids = get_background_bulbs_ids()

    print len(source_bulbs_ids)
    print len(background_bulbs_ids)
    # # building the interactome interface object
    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=False)
    local_matrix.full_rebuild()

    # # perform the interactome analysis
    # interactome_analysis([source_bulbs_ids], desired_depth=24, processors=3,
    #                      background_list=background_bulbs_ids, skip_sampling=False)

    # # building the reference parameters set
    # _filter = ['biological_process']
    # ref_param_set = [_filter, background_bulbs_ids, (1, 1), True, 3]

    # build the annotome interface
    # annot_matrix = AnnotomeInterface(*ref_param_set)
    # annot_matrix.full_rebuild()

    # perform the knowledge analysis
    # knowledge_analysis([source_bulbs_ids], desired_depth=24, processors=3,
    #                    param_set=ref_param_set, skip_sampling=False)
    # pass