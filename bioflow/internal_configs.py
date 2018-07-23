"""
Contains internal configurations that the user might want to act upon. For now, not documented in
the API_doc, just code comments for the maintainers
"""

#  Configures the mappings between concrete edge types and meta-types used for confidence
# calculation

## TODO: remove the edge_type_filter - it looks pretty deprecated by now.

edge_type_filters = {
    "Group": ["is_part_of_collection"],
    # Group relation group
    "Same": ["is_same"],
    # Same relation group
    "Reaction": ["is_Catalysant", "is_reaction_participant"],
    # Reaction relation group
    "Contact_interaction": ["is_part_of_complex", "is_Regulant"],
    # Contact_interaction relation group
    "HiNT_Contact_interaction": ["is_interacting"],
    # Contact_interaction relation group
    "BioGRID_Contact_interaction": ["is_weakly_interacting"],
    # possibly same interaction relations group
    "possibly_same": ["is_possibly_same"],
    # TF factors interaction
    "TRRUST_TF_Regulation": ["is_interacting"]}


#  Defines what nodes are to be masked to avoid conduction overload of non-informative nodes
Leg_ID_Filter = ['H+', 'ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2', 'NTP',
                 'Ubiquitin', 'cAMP', 'Actin']


#  Fudge for matrix diagonalization and other solver functions
line_loss = 1e-10

# Coefficients values for the value_Matrix
Adjacency_Martix_Dict = {
    "Group": 0.5,
    "Same": 1,
    "Reaction": 0.33,
    "Contact_interaction": 0.33,
    "weak_contact": 0.15,
    "possibly_same": 0.1, }

# Coefficients values for the conductance_Matrix
Conductance_Matrix_Dict = {
    "Group": 0.5,
    "Same": 100,
    "Reaction": 1,
    "Contact_interaction": 1,
    "weak_contact": 0.5,
    "possibly_same": 0.1, }

# allowed payload types for the annotation nodes
annotation_nodes_ptypes = [
    'name',
    'eCNumber',
    'ENSEMBL',
    'PathwayStep',
    'UniProt',
    'UNIPROT_Acnum',
    'UNIPROT_Name',
    'UNIPROT_GeneName',
    'UNIPROT_GeneRefs',
    'UNIPROT_GeneOL',
    'UNIPROT_GeneORF',
    'UNIPROT_Ensembl',
    'UNIPROT_EMBL_AC|~',
    'UNIPROT_EMBL_ID|~',
    'UNIPROT_PDB',
    'UNIPROT_GeneID', ]

interesting_lines = ['ID', 'AC', 'DE', 'GN', 'OX', 'DR']

interesting_xrefs = ['EMBL', 'GO', 'Pfam', 'Ensembl', 'KEGG', 'PDB', 'GeneID', 'SUPFAM']

names_to_ignore = [
    'Contains',
    'Allergen',
    'EC=',
    'Flags: ',
    'CD_antigen',
    'INN=']

uniprot_load_dict = {
    'Acnum': [],
    'Names': {
        'Full': '',
        'AltNames': []},
    'GeneRefs': {
        'Names': [],
        'OrderedLocusNames': [],
        'ORFNames': []},
    'Ensembl': [],
    'KEGG': [],
    'EMBL': [],
    'GO': [],
    'Pfam': [],
    'SUPFAM': [],
    'PDB': [],
    'GeneID': [],
    'RefSeq': [],
    'MGI': []}