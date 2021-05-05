"""
This is a set of top-level routines that have been wrapped for convenience
"""
from bioflow.annotation_network.BioKnowledgeInterface import \
    GeneOntologyInterface as AnnotomeInterface
from bioflow.molecular_network.InteractomeInterface import \
    InteractomeInterface as InteractomeInterface
from bioflow.neo4j_db.db_io_routines import cast_external_refs_to_internal_ids, \
    cast_background_set_to_bulbs_id, writer, Dumps
from bioflow.utils.io_routines import dump_object
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)


def map_and_save_gene_ids(hit_genes_location, all_detectable_genes_location=''):
    """
    Maps gene names/identifiers into internal database identifiers (neo4j ids) and saves them

    :param hit_genes_location: genes in the set we would like to analyse
    :param all_detectable_genes_location:  genes in the set that can be detected (background)
    :return: list of internal db ids for hits, list of internal db ids for background
    """

    standardized_hits = []  # [primary_set]
    standardized_secondary_hits = []  # [secondary_set=None]

    if type(hit_genes_location) == str:
        standardized_hits = [cast_external_refs_to_internal_ids(hit_genes_location)]
        standardized_secondary_hits = [None]

    if type(hit_genes_location) == tuple:
        standardized_hits = [cast_external_refs_to_internal_ids(hit_genes_location[0])]
        standardized_secondary_hits = [cast_external_refs_to_internal_ids(hit_genes_location[1])]

    if type(hit_genes_location) == list:
        for sub_hit_genes_location in hit_genes_location:
            if type(sub_hit_genes_location) == str:
                standardized_hits += [cast_external_refs_to_internal_ids(hit_genes_location)]
                standardized_secondary_hits += [None]
            if type(sub_hit_genes_location) == tuple:
                standardized_hits += [cast_external_refs_to_internal_ids(hit_genes_location[0])]
                standardized_secondary_hits += [cast_external_refs_to_internal_ids(hit_genes_location[1])]

    log.debug('standardized primary hits:\n\t%s' % standardized_hits)
    log.debug('standardized secondary_hits:\n\t%s' % standardized_secondary_hits)

    dump_object(Dumps.analysis_set_bulbs_ids, (standardized_hits, standardized_secondary_hits))

    if all_detectable_genes_location:
        background_set = cast_external_refs_to_internal_ids(all_detectable_genes_location)
        primary_set = [y for x in standardized_hits for y in x]
        formatted_secondary_hits = [_l
                                    if _l is not None
                                    else []
                                    for _l in standardized_secondary_hits]
        sec_set = [y for x in formatted_secondary_hits for y in x]

        primary_set = {_id
                       if (type(_id) == str, type(_id) == int)
                       else _id[0]
                       for _id in primary_set}

        sec_set = {_id
                   if (type(_id) == str, type(_id) == int)
                   else _id[0]
                   for _id in sec_set}

        if type(background_set[0]) == str or type(background_set[0]) == int:  # unweighted
            background_set = set(background_set).union(primary_set).union(sec_set)

        else:
            bck_set = {_id[0] for _id in sec_set}
            if not primary_set.issubset(bck_set):
                log.info('Nodes ids %s are missing in background set and are added with weight 0' %
                         (primary_set - bck_set))
                background_set += [(_id, 0) for _id in (primary_set - bck_set)]
            if not sec_set.issubset(bck_set):
                log.info('Secondary set nodes ids %s are missing in background set and are added '
                         'with weight 0' % (sec_set - bck_set))
                background_set += [(_id, 0) for _id in (sec_set - bck_set)]

    else:
        background_set = []

    dump_object(Dumps.background_set_bulbs_ids, background_set)

    # CURRENTPASS: save the background set as well as the active hits and secondary hits sets

    return standardized_hits, standardized_secondary_hits, background_set


def rebuild_the_laplacians():
    """
    Rebuilds the Annotome and Interactome interface objects in case of need,

    :return: None
    """
    local_matrix = InteractomeInterface()
    local_matrix.full_rebuild()

    annot_matrix = AnnotomeInterface()
    annot_matrix.full_rebuild()

