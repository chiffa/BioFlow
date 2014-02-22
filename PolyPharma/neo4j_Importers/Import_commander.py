__author__ = 'ank'
"""
Note: on a clean database, we should be able to load everything in ~ 6-7 hours.
Currently the bottleneck is on the UNIPROT indexation nodes
"""
####################################################################################
#
# Family of scripts regulating the whole import behavior
#
####################################################################################

from Reactome_org_inserter import clear_all, insert_all, run_diagnostics, full_dict
from GO_UNIPROT_Inserter import getGOs, import_GOs, import_UNIPROTS
from PolyPharma.neo4j_Declarations.General_operations import clean
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from Hint_importer import cross_ref_HiNT

# insert_all()
# run_diagnostics(full_dict)
# # clear_all(full_dict)
# import_GOs()
# # getGOs()
# # clean(DatabaseGraph.UNIPORT)
# import_UNIPROTS()
# # clean(DatabaseGraph.GOTerm)
# cross_ref_HiNT(True)

## run_diagnostics(full_dict)
