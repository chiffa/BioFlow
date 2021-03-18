"""
Contains internal configurations that the user might want to act upon.
"""

# REFACTOR: [better confs]: move to the configs .yaml
# Defines what nodes are to be masked to avoid conduction overload of non-informative nodes
# Their edges are not very hight necessarily, but they are highly central in the Reactome
# physical entities network
reactome_forbidden_nodes = ['H+', 'ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2',
                            'NTP', 'Ubiquitin', 'cAMP', 'Actin']

# Yeast uniprots that are excluded from calculation due to too high of edge values : 1300 - 3800
# Basically, members of protein degradation pathways
# CURRENTPASS: add human nodes to the neo4j parse
uniprot_forbidden_nodes = ['UBI4P_YEAST', 'NAB2_YEAST', 'SSB2_YEAST', 'HSP82_YEAST', 'SMT3_YEAST',
                           ]

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