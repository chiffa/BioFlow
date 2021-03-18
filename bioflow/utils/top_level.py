from bioflow.annotation_network.BioKnowledgeInterface import \
    GeneOntologyInterface as AnnotomeInterface
from bioflow.molecular_network.InteractomeInterface import \
    InteractomeInterface as InteractomeInterface
from bioflow.neo4j_db.db_io_routines import cast_analysis_set_to_bulbs_ids, \
    cast_background_set_to_bulbs_id, writer, Dumps
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids


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
        writer(open(Dumps.background_set_bulbs_ids, 'wt'), delimiter='\n').writerow(
            all_detectable_genes_ids)

    return hit_genes_ids, all_detectable_genes_ids


def rebuild_the_laplacians(all_detectable_genes=()):

    local_matrix = InteractomeInterface()
    local_matrix.full_rebuild()

    _filter = ['biological_process']
    # REFACTOR: [Better environment]: all_detectable_genes should not be here either
    ref_param_set = [_filter, all_detectable_genes, (1, 1), True, 3]

    annot_matrix = AnnotomeInterface(*ref_param_set)
    annot_matrix.full_rebuild()

