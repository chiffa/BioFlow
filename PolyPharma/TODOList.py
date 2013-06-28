'''
Created on Jun 25, 2013

@author: andrei
'''

# TODO: Unify DatabaseGraphDeclarations

# TODO: Include mappings to the Gene onthology
# TODO: Include mappings to the EC numbers from UNIPROT

# TODO: Refactor the UNIPROT parsing module
# It might be the most wise to first parse it with by specifying the tags you are interested in and then filling an Dictionary with them, then retrieving the tags one by one and storing 
# them in the wanted format (SQL table)
# On parse, add an additional "is_possibly the same" if a reactome protein displayName is the same as a UNIPROT SWISSPROT name without the 

# TODO: for each uniprot gene iserted add a DNA reference to the gene that was actually concerned.
# That sucks, because we will also need to insert all the names references and check for the possible collisions with the UNIPROT identifiers

# TODO: Matrix retriever:
# 1 -> 2 -> 3
#
# for each new insertion, look if it modifies any direct pathway
# even better: make all the insertions at 1 on the diagonal, then sup triangularise the matrix and fill in the same way as for the energy interaction: take the shortest interaction
# points (dynamic programming)
#
# Wait, can we do it by matrix multiplication? 
#
# 0 1 ?     0 1 2
#   0 1  ->   0 1
#     0         0

# We are not taking in account the EC numbers of enzymes, since they do not point towards the reaction participants, just the catalysing elements.