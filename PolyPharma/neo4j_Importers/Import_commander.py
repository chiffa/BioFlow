__author__ = 'ank'

####################################################################################
#
# Family of scripts regulating the whole import behavior
#
####################################################################################

from Reactome_org_inserter import clear_all, insert_all, run_diagnostics, full_dict

#insert_all()
run_diagnostics(full_dict)
#clear_all(full_dict)

