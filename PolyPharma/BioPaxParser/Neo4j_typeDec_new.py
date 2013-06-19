'''
Created on Jun 13, 2013

@author: andrei
'''

from bulbs.model import Node, Relationship
from bulbs.property import String, Integer, Float, Dictionary, List
from bulbs.utils import current_datetime

class CostumNode(Node):             # Serves as a basis for the annotation
    element_type="CostumNode"
    ID=String(nullable=False)       # Reactome import heritage
    displayName=String()            # To see what it is, for the human operator
    custom=String()                 # Just in case
    load=Float()                    # To freeze information transmission score (Dict should be better?)

class AnnotNode(Node):                  # Used mainly the simplest annotation basis annotation
    element_type="AnnotNode"   
    ContentType=String(nullable=False)  # Payload type          
    Contents=String(nullable=False)     # Indexed payload
 
class CostumLink(Relationship):
    label="CostumLink"
    linkType=String()
    custom=String()
    load=Float()

#<====================================>
    
class Meta(CostumNode):
    element_type="Meta"
    
class Fragment(CostumNode):
    element_type="Fragment"
    locationType=String()           # SPAN or POINT
    location=String(nullable=False) # Location is protein-relative

class Instantiation_Type(CostumNode):
    element_type="Instantiation_Type"
    type=String()                   # Instantiation type: 

class Instance(CostumNode):
    element_type="Instance"
    
class Collection(CostumNode):
    element_type="Collection"

class Reaction(CostumNode):
    element_type="Reaction"
    frequency=Float()
    
#<=====================================>

class DNA(Meta):
    element_type="DNA"

class Location(Instantiation_Type):
    element_type="Location"

class is_localized(Instantiation_Type):
    element_type="is_localized"

class DNA_Collection(Meta):
    element_type="DNA Collection"

class is_part_of_collection(CostumLink):
    label="is_part_of_collection"

class is_annotated(Relationship):
    label="is_annotated"

class Complex(Meta):
    element_type="Complex"

class is_part_of_complex(Relationship):
    label="is_part_of_complex"

class Complex_Collection(Meta):
    element_type="Complex_Collection"

class is_catalysant(CostumLink):
    label="is_Catalysant"
    controlType=String()
    ID=String(nullable=False)
    displayName=String()
    
class is_regulant(CostumLink):
    label="is_Regulant"
    controlType=String()
    ID=String(nullable=False)
    displayName=String()

class PhysicalEntity(Meta):
    element_type="PhysicalEntity"
    
class PhysicalEntity_Collection(Meta):
    element_type="PhysicalEntity_Collection"

class TemplateReaction(Reaction):
    element_type="TemplateReaction"
    
class Degradation(Reaction):
    element_type="Degradation"

class RNA(Meta):
    element_type="RNA"
    
class RNA_Collection(Meta):
    element_type="RNA_Collection"

class Originating_Organism(Instantiation_Type):
    element_type="Originating_Organism"

class is_originating_in_organism(Relationship):
    label="is_originating_in_organism"
    
class Protein(Meta):    
    element_type="Protein"

class Protein_Collection(Meta):
    element_type="Protein_Collection"

class SmallMolecule(Meta):
    element_type="Small_Molecule"

class SmallMolecule_Collection(Meta):
    element_type="Small_Molecule_Collection"

class BiochemicalReaction(Reaction):
    element_type="Biochemical_Reaction"
    
class is_reaction_participant(CostumLink):
    label="is_reaction_particpant"
    side=String()

class ModificationFeature(Instantiation_Type):
    element_type="Modification_Feature"
    location=String()

class is_modified_to(CostumLink):
    label="is_modified_to "

class is_able_to_modify(CostumLink):
    label="is_able_to_modify"


        