"""
This module contains all the declaration of object models representing the underlying database state
Currently all the declaration is done in Python, which leads to a plethora of shallow classes
disliked by PyLint, but currently preferred to declaration through a meta-language.
"""

# WARNING FOR MAINTAINERS: There should be no difference in the "element_type" or "label" attributes
# compared to the class name. They are used in the indexes and metaprogramming specifying how the
# names of nodes map to neo4j's node indexes

from bulbs.model import Node, Relationship
from bulbs.property import String, Integer, Float, Bool


class CostumNode(Node):             # Serves as a basis for the annotation
    element_type = "CostumNode"
    ID = String(nullable=False)
    # TODO: rename to the legacy ID. This is a Reactome import legacy
    displayName = String()            # To see what it is, for the human operator
    main_connex = Bool()
    custom = String()                 # Just in case
    load = Float()                    # Deprecated. To freeze information transmission score
    # TODO: remove deprecated properties and rename properties their wanted
    # names


class AnnotNode(Node):                # Used mainly the simplest annotation basis annotation
    """
    Available AnnotNode Types are listed in AnnotNode_ptypes
    """
    element_type = "AnnotNode"
    ptype = String(nullable=False)    # Payload type
    payload = String(nullable=False)  # Indexed payload
    # TODO: move the supported annotation types nodes to the main
    # documentation set.


class CostumLink(Relationship):
    label = "CostumLink"
    costum_from = String()
    costum_to = String()
    linkType = String()
    custom = String()
    load = Float()
    source = String()
    weight = Float()
# < ======================================================================= >


class Meta(CostumNode):
    element_type = "Meta"
    localization = String()


class Fragment(CostumNode):
    element_type = "Fragment"
    locationType = String()             # SPAN or POINT
    location = String(nullable=False)   # Location is protein-relative


class Instantiation_Type(CostumNode):
    element_type = "Instantiation_Type"
    type = String()                   # Instantiation type:


class Instance(CostumNode):
    element_type = "Instance"


class Collection(CostumNode):
    element_type = "Collection"


class Reaction(CostumNode):
    element_type = "Reaction"
    frequency = Float()
# < ========================================================================= >


class DNA(Meta):
    element_type = "DNA"


class Location(Instantiation_Type):
    element_type = "Location"


class is_localized(CostumLink):
    label = "is_localized"


class DNA_Collection(Meta):
    element_type = "DNA_Collection"


class is_part_of_collection(CostumLink):
    label = "is_part_of_collection"


class is_annotated(CostumLink):
    label = "is_annotated"


class is_possibly_same(CostumLink):
    label = "is_possibly_same"


class Complex(Meta):
    element_type = "Complex"


class is_part_of_complex(CostumLink):
    label = "is_part_of_complex"


class Complex_Collection(Meta):
    element_type = "Complex_Collection"


class is_catalysant(CostumLink):
    label = "is_catalysant"  # correct here!
    controlType = String()
    ID = String(nullable=False)
    displayName = String()


class is_regulant(CostumLink):
    label = "is_regulant"  # correct here!
    controlType = String()
    ID = String(nullable=False)
    displayName = String()


class PhysicalEntity(Meta):
    element_type = "PhysicalEntity"


class PhysicalEntity_Collection(Meta):
    element_type = "PhysicalEntity_Collection"


class TemplateReaction(Reaction):
    element_type = "Template_Reaction"


class Degradation(Reaction):
    element_type = "Degradation"


class RNA(Meta):
    element_type = "RNA"


class RNA_Collection(Meta):
    element_type = "RNA_Collection"


class Originating_Organism(Instantiation_Type):
    element_type = "Originating_Organism"


class is_originating_in_organism(Relationship):
    label = "is_originating_in_organism"


class Protein(Meta):
    element_type = "Protein"


class Protein_Collection(Meta):
    element_type = "Protein_Collection"


class SmallMolecule(Meta):
    element_type = "SmallMolecule"  # correct here! (?)


class SmallMolecule_Collection(Meta):
    element_type = "SmallMolecule_Collection"  # correct here!


class BiochemicalReaction(Reaction):
    element_type = "BiochemicalReaction"  # correct here!


class is_reaction_participant(CostumLink):
    label = "is_reaction_particpant"
    side = String()


class ModificationFeature(Instantiation_Type):
    element_type = "ModificationFeature"  # correct here!
    location = String()


class is_modified_to(CostumLink):
    label = "is_modified_to "


class is_modified_by(CostumLink):
    label = "is_modified_by "


class is_able_to_modify(CostumLink):
    label = "is_able_to_modify"


class Pathway(CostumNode):
    element_type = "Pathway"


class Pathway_Step(CostumNode):
    element_type = "Pathway_Step"


class is_part_of_pathway(CostumLink):
    label = "is_part_of_pathway"


class is_next_in_pathway(CostumLink):
    label = "is_next_in_pathway"


class UNIPROT(Meta):
    element_type = "UNIPROT"
    # 1 / 0 -> if it is involved in the Reactome.org interactions or not.
    involved = Bool()


class GOTerm(CostumNode):
    element_type = "GOTerm"
    Name = String()
    Namespace = String()
    Definition = String()


class is_go_annotation(CostumLink):
    label = "is_go_annotation"


class is_a_go(CostumLink):
    label = "is_a_go"


class is_part_of_go(CostumLink):
    label = "is_part_of_go"


class is_same(CostumLink):
    label = "is_same"


class is_interacting(
        CostumLink):   # According to the Yu lab hint database for humans there is one-to-one interaction
    label = "is_interacting"


class is_weakly_interacting(CostumLink):  # interaction according to biogrid_path
    label = "is_weakly_interacting"
    throughput = String()      # high or all
    confidence = Float()        # float
