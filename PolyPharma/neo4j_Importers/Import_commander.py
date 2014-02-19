__author__ = 'ank'

####################################################################################
#
# Family of scripts regulating the whole import behavior
#
####################################################################################

from Reactome_org_inserter import clear_all, insert_all, run_diagnostics, full_dict
from GO_UNIPROT_Inserter import clean, getGOs, import_GOs, import_UNIPROTS
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from Hint_importer import cross_ref_HiNT

# insert_all()
# run_diagnostics(full_dict)
# # clear_all(full_dict)
# import_GOs()
# # getGOs(DatabaseGraph.GOTerm)
# import_UNIPROTS()
# # clean(DatabaseGraph.UNIPORT)
# # clean(DatabaseGraph.GOTerm)
cross_ref_HiNT(True)