from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph

forbidden_verification_dict = {
    'Small Molecule': DatabaseGraph.SmallMolecule,
    'Small Molecule Collection': DatabaseGraph.SmallMolecule_Collection,
    'Physical Entity': DatabaseGraph.PhysicalEntity,
    'Physical Entity Collection': DatabaseGraph.PhysicalEntity_Collection,
    }

full_dict = {
    'DNA': (DatabaseGraph.DNA, "DNA"),
    'DNA Collection': (DatabaseGraph.DNA_Collection, "DNA_Collection"),
    'RNA': (DatabaseGraph.RNA, "RNA"),
    'RNA Collection': (DatabaseGraph.RNA_Collection, "RNA_Collection"),
    'Small Molecule': (DatabaseGraph.SmallMolecule, "SmallMolecule"),
    'Small Molecule Collection': (DatabaseGraph.SmallMolecule_Collection, "SmallMolecule_Collection"),
    'Protein': (DatabaseGraph.Protein, "Protein"),
    'Protein Collection': (DatabaseGraph.Protein_Collection, "Protein_Collection"),
    'Complex': (DatabaseGraph.Complex, "Complex"),
    'Complex Collection': (DatabaseGraph.Complex_Collection, "Complex_Collection"),
    'Physical Entity': (DatabaseGraph.PhysicalEntity, "PhysicalEntity"),
    'Physical Entity Collection': (DatabaseGraph.PhysicalEntity_Collection, "PhysicalEntity_Collection"),
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
    }
