"""
Contains internal configurations that the user might want to act upon. For now, not documented in
the API_doc, just code comments for the maintainers
"""

# CURRENTPASS: [EDGE WEIGHTS] refactor this into a more clear state.

#  Configures the mappings between concrete edge types and meta-types used for confidence
# calculation

# CURRENTPASS: Deprecate this
# => replace the
deprecated_edge_type_filters = {
    "Group": ["is_part_of_collection"],
    # Group relation group
    "Same": ["is_same"],
    # Same relation group
    "Reaction": ["is_catalysant", "is_reaction_participant"],
    # Reaction relation group
    "Contact_interaction": ["is_part_of_complex", "is_regulant"],
    # Contact_interaction relation group
    "HiNT_Contact_interaction": ["is_interacting"],
    # Contact_interaction relation group
    "BioGRID_Contact_interaction": ["is_weakly_interacting"],
    # possibly same interaction relations group
    "possibly_same": ["is_likely_same"],
    # TF factors interaction
    "TRRUST_TF_Regulation": ["is_interacting"]}


# CURRENTPASS: move to the configs .yaml
# Defines what nodes are to be masked to avoid conduction overload of non-informative nodes
# Their edges are not very hight necessarily, but they are highly central in the Reactome
# physical entities network
reactome_forbidden_nodes = ['H+', 'ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2',
                            'NTP', 'Ubiquitin', 'cAMP', 'Actin']

# Yeast uniprots that are excluded from calculation due to too high of edge values : 1300 - 3800
# Basically, members of protein degradation pathways
# TODO: add human nodes to the neo4j parse
uniprot_forbidden_nodes = ['UBI4P_YEAST', 'NAB2_YEAST', 'SSB2_YEAST', 'HSP82_YEAST', 'SMT3_YEAST']


# Coefficients values for the value_Matrix
deprecated_adjacency_matrix_weights = {
    "Group": 0.5,
    "Same": 1,
    "Reaction": 0.33,
    "Contact_interaction": 0.33,
    "weak_contact": 0.15,
    "is_likely_same": 0.1,
    "skipped": 0}

# Coefficients values for the conductance_Matrix
deprecated_laplacian_matrix_weights = {
    "Group": 0.5,
    "Same": 100,
    "Reaction": 1,
    "Contact_interaction": 1,
    "weak_contact": 0.5,
    "is_likely_same": 1,
    "skipped": 0}

# allowed payload types for the annotation nodes
# CURRENTPASS: not used anymore
deprecated_annotation_nodes_ptypes = [
    'name',
    'eCNumber',
    'ENSEMBL',
    'PathwayStep',
    'UniProt',
    'UNIPROT_Acnum',
    'UNIPROT_Name',
    'UNIPROT_AltName'
    'UNIPROT_GeneName',
    'UNIPROT_AltGeneName',
    'UNIPROT_GeneRefs',
    'UNIPROT_GeneOL',
    'UNIPROT_GeneORF',
    'UNIPROT_Ensembl',
    'UNIPROT_EMBL_AC|~',
    'UNIPROT_EMBL_ID|~',
    'UNIPROT_PDB',
    'UNIPROT_GeneID', ]

# CURRENTPASS: Investigate the use
#  Eithere a dynamic registration (in which case Reactome importing is unskippable)
#  > Or a recovery from the neo4j as "source:Reactome"
#  Or user configuration file of what to import
to_deprecate_neo4j_names_dict = {
    'DNA': "DNA",
    'DNA Collection': "DNA_Collection",
    'RNA': "RNA",
    'RNA Collection': "RNA_Collection",
    'Small Molecule': "SmallMolecule",
    'Small Molecule Collection': "SmallMolecule_Collection",
    'Protein': "Protein",
    'Protein Collection': "Protein_Collection",
    'Complex': "Complex",
    'Complex Collection': "Complex_Collection",
    'Physical Entity': "PhysicalEntity",
    'Physical Entity Collection': "PhysicalEntity_Collection",
    'TemplateReaction': "Template_Reaction",
    'Degradation': "Degradation",
    'BiochemicalReaction': "BiochemicalReaction",
    'Pathway Step': "PathwayStep",
    'Pathway': "Pathway",
    'Cell Locations': "Location",
    'Annotations': "AnnotNode",
    'Modification Feature': "ModificationFeature",
    'UNIPROT': "UNIPROT",
    'GO Term': "GOTerm",
    'COMPLEX': "COMPLEX",
}

full_list = list(to_deprecate_neo4j_names_dict.keys())

# CURRENTPASS: deprecate them
deprecated_reactome_reactions_types_list = ['TemplateReaction', 'Degradation', 'BiochemicalReaction']