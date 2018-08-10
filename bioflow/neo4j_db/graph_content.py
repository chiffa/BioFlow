
neo4j_names_dict = {
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
    'Pathway Step': "Pathway_Step",
    'Pathway': "Pathway",
    'Cell Locations': "Location",
    'Annotations': "AnnotNode",
    'Modification Feature': "ModificationFeature",
    'UNIPROT': "UNIPROT",
    'GO Term': "GOTerm",
    'COMPLEX': "COMPLEX",
}


forbidden_verification_list = ['Small Molecule',
                               'Small Molecule Collection',
                               'Physical Entity',
                               'Physical Entity Collection',
                               'Protein',
                               'Protein Collection',
                               'UNIPROT'
                               ]


full_list = list(neo4j_names_dict.keys())
