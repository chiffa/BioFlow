"""
Contains internal configurations that the user might want to act upon. For now, not documented in
the API_doc, just code comments for the maintainers
"""

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

# # allowed payload types for the annotation nodes
# deprecated_annotation_nodes_ptypes = [
#     'name',
#     'eCNumber',
#     'ENSEMBL',
#     'PathwayStep',
#     'UniProt',
#     'UNIPROT_Acnum',
#     'UNIPROT_Name',
#     'UNIPROT_AltName'
#     'UNIPROT_GeneName',
#     'UNIPROT_AltGeneName',
#     'UNIPROT_GeneRefs',
#     'UNIPROT_GeneOL',
#     'UNIPROT_GeneORF',
#     'UNIPROT_Ensembl',
#     'UNIPROT_EMBL_AC|~',
#     'UNIPROT_EMBL_ID|~',
#     'UNIPROT_PDB',
#     'UNIPROT_GeneID', ]

# CURRENTPASS: Investigate the use
#  Eithere a dynamic registration (in which case Reactome importing is unskippable)
#  > Or a recovery from the neo4j as "source:Reactome"
#  Or user configuration file of what to import
to_deprecate_neo4j_names_dict = {  # TRACING: [deprecation]
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
    'UNIPROT': "GOTerm",
    'GO Term': "GOTerm",
    'COMPLEX': "COMPLEX",
}