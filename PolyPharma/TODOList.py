'''
Created on Jun 25, 2013

@author: andrei
'''

# TODO: reimplement the whole GO program based on the matrixes and the correct information circulation attributes
# TODO: do clustering by the method suggested by Misiuro: we might want clusters that contain a particular target
    # several times, since it might contribute to several different targets


# Cancelled: Include mappings to the EC numbers from UNIPROT

# TODO: possible extensions to the organ specificity

# TODO: (See with Bourne): perform a HiNT-only analysis of reactivity

# Cancelled (done on pulling): add a filter on the GO import to exclude:
    # - Molecular functions
    # - cellular locations

# TODO: add pdb references to all the indexed proteins from uniprot

# TODO: for each uniprot gene iserted add a DNA reference to the gene that was actually concerned.
# That sucks, because we will also need to insert all the names references and check for the possible collisions with the UNIPROT identifiers


# We are not taking in account the EC numbers of enzymes, since they do not point towards the reaction participants, just the catalysing elements.

# TODO: Ask Li if the eHiTs scores are proportional to the binding energy / inhibition strength

# Ok, so the method is assymetric, which is really uncool

# TODO: once set and intermediate term construction is terminated, we would need to
        # create a new dictionary containing the GO terms to informativity mapping
        # The informativity transformation should be done on loading. 
        # In fact this informativity transformation would allow us to balance the 
        # term informativity (what is the fraction of proteins that have it)
        # v.s. the number that this GO term appears among different terms of a 
        # group of terms
        
# TODO: Implement a benchmark comparing the global term informativity v.s. the local
        # informativity and see how a flow methods relates to the benchmark test.
        
        ###
        # Well, what happens is that we are interested in computing the information
        # flow, and not the gradient that induces that flow. So our GO are resistances
        # but our proteins are in fact sources of current and not of tension
        
        ###
        # The problem is the management of the forced current flows. It is much better
        # to consider it as a system with total tension proportional to the importance
        # of all the proteins loaded and a conductivity of each protein and GO term
        # that is proportional to it's importance or informativeness
        
        ###
        # An alternative would be to use Lagrange method for optimization under
        # Under Constraints, where the constraint is actaully the information
        # flow through each of the proteins
        
        ####
        # Wouldn't it be easier to use java-like object oriented programming instead
        # of using a set of inter-indexed dictionnaries?
        # The problem is that the dictionaries are not necessarily visible
        # On the other hand, we might have to manage a non-shallow copies of 
        # objects in case they are 
        # Err. I guess after all it will actually be easier to implement it all with
        # Proteins as objects and GOs as objects, which are created individually
        # for each class.

        ####
        # TODO: build a benchmark to compare the results for the two-sided and one
        # sided protein computation. My guess is that the two-sided would be 
        # slightly more sensitive to the weights of proteins
        
        ####
        # TODO: build a benchmark based on LU transformation and information flow based
        # on a small subset of proteins. Each protein is a Node, each GO is an edge

# Ideas regarding the information analysis: get the most information circulation within group and the most different terms for inter-group information circulation

