"""
Contains internal configurations that the user might want to act upon. For now, not documented
"""
__author__ = 'ank'

########################################################################################################
#  Configures the mappings between concrete edge types and meta-types used for confidence computation  #
########################################################################################################
edge_type_filters = {
    "Group" : ["is_part_of_collection"],                                  # Group relation group
    "Same" : ["is_same"],                                                 # Same relation group
    "Reaction" : ["is_Catalysant", "is_reaction_particpant"],             # Reaction relation group
    "Contact_interaction" : ["is_part_of_complex", "is_Regulant"],        # Contact_interaction relation group
    "HiNT_Contact_interaction" : ["is_interacting"],                      # Contact_interaction relation group
    "BioGRID_Contact_interaction": ["is_weakly_interacting"],
    "possibly_same" : ["is_possibly_same"],
    }

##############################################################################################
#  Defines what nodes are to be masked to avoid conduction overload of non-informative nodes #
##############################################################################################
Leg_ID_Filter = ['H+', 'ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2', 'NTP',]

############################################################################
#  Fundge for matrix diagolizations of matrixes and other solver functions #
############################################################################
fudge = 1e-10

# Coefficients values for the value_Matrix
Adjacency_Martix_Dict = {"Group":0.5,
             "Same":1,
             "Reaction":0.33,
             "Contact_interaction":0.33,
             "weak_contact": 0.15,
             "possibly_same":0.1,
             }


# Coefficients values for the conductance_Matrix
Conductance_Matrix_Dict = {"Group":0.5,
             "Same":100,
             "Reaction":1,
             "Contact_interaction":1,
             "weak_contact":0.5,
             "possibly_same":0.1,
             }

