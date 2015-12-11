"""
The main method of the neo4j database building

On a clean database, we should be able to yeast databases in ~ 6-7 hours and humans in less than
24 hours.
"""
from BioFlow.neo4j_db.db_io_routines import recompute_forbidden_IDs, Forbidden_verification_dict
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph
from BioFlow.db_importers.go_and_uniprot_importer import getGOs, import_GOs, import_UNIPROTS
from BioFlow.db_importers.hint_importer import cross_ref_hint
from BioFlow.db_importers.reactome_importer import clear_all, insert_all, run_diagnostics, full_dict
from biogrid_importer import import_bio_grid


# TODO: add the aboundance import
# TODO: add the derivative importance contribution

# #################################
# redirecting all to a log file
# f = open('../logs/Commander_logs.log','w')
# sys.stdout = f
# ################################

def build_db():
    insert_all()
    import_GOs()
    import_UNIPROTS()
    cross_ref_hint(flush=True)
    import_bio_grid()
    run_diagnostics(full_dict)
    recompute_forbidden_IDs(Forbidden_verification_dict)


def destroy_db():
    clear_all(full_dict)


if __name__ == "__main__":
    # clear_all(full_dict)
    # run_diagnostics(full_dict)
    # insert_all()
    # run_diagnostics(full_dict)

    # clean(DatabaseGraph.GOTerm)
    # import_GOs()

    getGOs()
    clear_all({'UNIPROT': (DatabaseGraph.UNIPORT, "UNIPROT")})
    import_UNIPROTS()

    cross_ref_hint()

    import_bio_grid()

    run_diagnostics(full_dict)

    recompute_forbidden_IDs(Forbidden_verification_dict)

    pass
