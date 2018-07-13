from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph, on_alternative_graph

if on_alternative_graph:
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

else:
    neo4j_names_dict = {
        'DNA': (DatabaseGraph.DNA, "DNA"),
        'DNA Collection': (DatabaseGraph.DNA_Collection, "DNA_Collection"),
        'RNA': (DatabaseGraph.RNA, "RNA"),
        'RNA Collection': (DatabaseGraph.RNA_Collection, "RNA_Collection"),
        'Small Molecule': (DatabaseGraph.SmallMolecule, "SmallMolecule"),
        'Small Molecule Collection': (DatabaseGraph.SmallMolecule_Collection,
                                      "SmallMolecule_Collection"),
        'Protein': (DatabaseGraph.Protein, "Protein"),
        'Protein Collection': (DatabaseGraph.Protein_Collection, "Protein_Collection"),
        'Complex': (DatabaseGraph.Complex, "Complex"),
        'Complex Collection': (DatabaseGraph.Complex_Collection, "Complex_Collection"),
        'Physical Entity': (DatabaseGraph.PhysicalEntity, "PhysicalEntity"),
        'Physical Entity Collection': (DatabaseGraph.PhysicalEntity_Collection,
                                       "PhysicalEntity_Collection"),
        'TemplateReaction': (DatabaseGraph.TemplateReaction, "Template_Reaction"),
        'Degradation': (DatabaseGraph.Degradation, "Degradation"),
        'BiochemicalReaction': (DatabaseGraph.BiochemicalReaction, "BiochemicalReaction"),
        'Pathway Step': (DatabaseGraph.PathwayStep, "Pathway_Step"),
        'Pathway': (DatabaseGraph.Pathway, "Pathway"),
        'Cell Locations': (DatabaseGraph.Location, "Location"),
        'Annotations': (DatabaseGraph.AnnotNode, "AnnotNode"),
        'Modification Feature': (DatabaseGraph.ModificationFeature, "ModificationFeature"),
        'UNIPROT': (DatabaseGraph.UNIPORT, "UNIPROT"),
        'GO Term': (DatabaseGraph.GOTerm, "GOTerm"),
        'COMPLEX': (DatabaseGraph.COMPLEX, "COMPLEX"),
        }


# TODO: perform an edit here to kill the UBC insanity
forbidden_verification_list = ['Small Molecule',
                               'Small Molecule Collection',
                               'Physical Entity',
                               'Physical Entity Collection',
                               'Protein',
                               'Protein Collection',
                               'UNIPROT'
                               ]

full_list = list(neo4j_names_dict.keys())
