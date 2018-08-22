"""
The main method of the neo4j database building

On a clean database, we should be able to yeast databases in ~ 6-7 hours and humans in less than
24 hours.
"""
from bioflow import main_configs
from bioflow.bio_db_parsers.geneOntologyParser import GOTermsParser
from bioflow.bio_db_parsers.uniprotParser import UniProtParser
from bioflow.db_importers.hint_importer import cross_ref_hint
from bioflow.db_importers.reactome_importer import insert_reactome
from bioflow.db_importers.biogrid_importer import cross_ref_bio_grid
from bioflow.db_importers.tf_importers import cross_ref_tf_factors
from bioflow.db_importers.phosphosite_importer import cross_ref_kinases_factors
from bioflow.db_importers.complex_importer import insert_complexes
from bioflow.db_importers.go_and_uniprot_importer import memoize_go_terms, import_gene_ontology, \
    import_uniprots, pull_up_acc_nums_from_reactome
from bioflow.neo4j_db.db_io_routines import recompute_forbidden_ids, clear_all, run_diagnostics,\
    cross_link_identifiers, compute_annotation_informativity
from bioflow.neo4j_db.graph_content import forbidden_verification_list, full_list


def build_db():
    insert_reactome()

    _go_terms, _go_terms_structure = GOTermsParser().parse_go_terms(main_configs.gene_ontology_path)
    import_gene_ontology(_go_terms, _go_terms_structure)

    _uniprot = UniProtParser(main_configs.up_tax_ids).parse_uniprot(main_configs.uniprot_path)
    _reactome_acnum_bindings = pull_up_acc_nums_from_reactome()
    import_uniprots(_uniprot, _reactome_acnum_bindings)

    cross_ref_hint()
    cross_ref_bio_grid()

    cross_ref_tf_factors('t')
    cross_ref_kinases_factors()
    insert_complexes()

    run_diagnostics(full_list)

    cross_link_identifiers()

    compute_annotation_informativity()

    recompute_forbidden_ids(forbidden_verification_list)


def destroy_db():
    clear_all(full_list)


if __name__ == "__main__":
    pass
    # clear_all(full_list)
    # run_diagnostics(full_list)
    # insert_reactome()
    # run_diagnostics(full_list)
    #
    # clear_all(['GO Term'])
    #
    # go_terms, go_terms_structure = GOTermsParser().parse_go_terms(main_configs.gene_ontology_path)
    # import_gene_ontology(go_terms, go_terms_structure)

    # memoize_go_terms()
    #
    # clear_all(['UNIPROT'])
    #
    # uniprot = UniProtParser(main_configs.up_tax_ids).parse_uniprot(main_configs.uniprot_path)
    #
    # reactome_acnum_bindings = pull_up_acc_nums_from_reactome()
    # import_uniprots(uniprot, reactome_acnum_bindings)
    #
    # cross_ref_hint()
    # cross_ref_bio_grid()
    #
    # cross_ref_tf_factors('t')
    #
    # cross_ref_kinases_factors()
    #
    # insert_complexes()
    #
    # run_diagnostics(full_list)
    #
    # cross_link_identifiers()
    #
    # compute_annotation_informativity()
    #
    recompute_forbidden_ids(forbidden_verification_list)
