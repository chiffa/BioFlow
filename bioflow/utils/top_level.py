"""
This is a set of top-level routines that have been wrapped for convenience
"""
from bioflow.annotation_network.BioKnowledgeInterface import \
    GeneOntologyInterface as AnnotomeInterface
from bioflow.molecular_network.InteractomeInterface import \
    InteractomeInterface as InteractomeInterface
from bioflow.neo4j_db.db_io_routines import cast_analysis_set_to_bulbs_ids, \
    cast_background_set_to_bulbs_id, writer, Dumps
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids


def map_and_save_gene_ids(hit_genes_location, all_detectable_genes_location=''):
    """
    Maps gene names/identifiers into internal database identifiers (neo4j ids) and saves them

    :param hit_genes_location: genes in the set we would like to analyse
    :param all_detectable_genes_location:  genes in the set that can be detected (background)
    :return: list of internal db ids for hits, list of internal db ids for background
    """
    cast_analysis_set_to_bulbs_ids(hit_genes_location)
    hit_genes_ids = get_source_bulbs_ids()
    print('debug, top_level hit_genes_ids: %s' % hit_genes_ids)

    if all_detectable_genes_location:
        cast_background_set_to_bulbs_id(
            background_set_csv_location=all_detectable_genes_location,
            analysis_set_csv_location=hit_genes_location)

        all_detectable_genes_ids = get_background_bulbs_ids()

    else:
        all_detectable_genes_ids = []
        writer(open(Dumps.background_set_bulbs_ids, 'wt'), delimiter='\n').writerow(
            all_detectable_genes_ids)

    return hit_genes_ids, all_detectable_genes_ids


def rebuild_the_laplacians():
    """
    Rebuilds the Annotome and Interactome interface objects in case of need,

    :return: None
    """
    local_matrix = InteractomeInterface()
    local_matrix.full_rebuild()

    annot_matrix = AnnotomeInterface()
    annot_matrix.full_rebuild()

