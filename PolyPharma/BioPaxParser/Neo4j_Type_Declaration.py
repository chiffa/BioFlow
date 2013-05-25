'''
Created on 13 mai 2013

@author: Andrei Kucharavy

TODO: learn how to build indexes with bulbs
TODO: see what is accessed when a index is an exact hit. iterator? List?

 Can we do class heritages? Like the ones in the owl system: localized version comprehends the same
 properties as the non-localized one, for instance modification features, domains, conserved regions,
 splicing variants, instances included in compounds.
 
 1) Create separate instanciated features => potentially exponential node growth, at the limit of neo4j 
 capacities
 2) Information expressed in the links that are intercating: problem: information extraction about a target
 and predictive work might get VEEERY complicated
 3) Lazy instance creation: create only if scientific evidence already exists or if the given features are
 critical for the interaction: originally one instanciated pull of molecules
    =?> how to conserve all the possible modifications?
       =!> Create an instance for each possible distinct instanciation event and then point the possible instantiation
       events from the instances critical for the interaction
            =!> instantiation_type; properties: permitted/forbidden,  inst_type, inst_location
            =!> instantiation
 
 
 In fact I am trying to create a sufficiently powerful machine-readable language allowing to describe the properties
 of interesting biological molecules to a single molecule level.
 A coding language for biology able to describe anything said in current research
 
 In a more applied fashion: does UNIPROT distinguish proteins issued from alternative splicing forms? Any other database?
 
 Proof that a possible alternative splicing actually exists. Fuck, splicing isn't even in the BioPax definition. This sucks a lot
 
 Catalytic site activity modulation by three most common post-translational modifications
 
 The most promising approach is to start off with genes, then define from them differently spliced translation and transcription;
 Then define mutated versions of genes regarding the 
 
 for the meta-regulation relations: three sub-types
     - regulate/modify
     - can regulate/modify
     - Cannot regulate/modify
     
 Since we are interested mainly in the information transition, any specific relations are right now not necessary
 
 * "Small Molecule" is just anything which is a DNA, RNA, Protein, or Complex
 
 Critical processes/compounds: Meta-Level
 
     - DNA packaging / unpackaging
     - RNA transcription
     - Protein translation
     - Compex assembly
     
     - Protein catalytic activity
     - RNA catalytic activity
     - Complex catalytic activity
     - "Small Molecule" catalytic activity
     
     - Protein processing (post-translational modification, expoort, import, degradation, conformation modification)
     - Small molecule processing
     - RNA processing
     - DNA processing (includes spontaneous modification, retro-transposon insertion, etc...)
     - Complex processing (post-translational modification, )
     
     - Reactions:
        - Transport
        - post-translational modification
        - Conformation modification (Instantiation modification)
        - Chemical reaction sensu strictu (processing of small molecules/ transforming proteins one into an other)
        - Protein degradation
        - DNA/RNA Alteration
          -> Insertion 
          -> Excision
          -> Point-mutation
        - DNA packaging modification
        - mRNA alternative splicing
        - Degradation
     
     # Fragmentation types
     
     - DNA fragmentation:
        - Gene
          - exon/intron
        - Regulatory element
        - Conserved domain  \ unite them together by definig an "evolution speed-defined element"
        - evolution hotspot /
        - Retrotransposon / other element of interest
        - Protein binding site
        - RNA binding site
        - encoded conserved region
        - encoded part of catalytic region
        - encoded post-translation modification site
        - encoded domain
    - RNA fragmentation
        - exon/intron
        - catalytic region / protein-non-relevant domain
        - conserved domain  \ unite them together by definig an "evolution speed-defined element"
        - evolution hotspot /
        - protein binding region
        - complex binding region
        - small molecule binding region
        - encoded conserved region
        - encoded part of catalytic region
        - encoded post-translation modification site
        - encoded domain
    - Protein fragmentation
        - conserved domain  \ unite them together by definig an "evolution speed-defined element" => include swipes, 
        - evolution hotspot /
        - part of catalytic site
        - post-translational modification site
        - protein domain
        - catalytic site
    - Complex fragmentation:
        - Proteins
        - RNA
        - Small molecule
        - Composite catalytic sites (i.e. catalytic sites that doesn't exist outside of the complex)
        
    
    # Instantiation event types: 
    - Localization in a cellular compartment
    - Mutation
    - Post-translational modification 
    - NOT alternative splicing: treated as transcription to different RNA
    - Participation in a complex
    - Conformation modification
    
    
 
 
 
 Major Model Repositry and associated search method:
    - resolution of the query terms (typing and accessing)
    - resolution of the query term correlation
    - mapping of the correlation as the shortest term 
    - regulation of depths of search
 
 degradation might also be an interesting action to add. 
 
 TODO: add an evidence describing method for the observation of different events?
 TODO: catalytic action of messenger RNA?
 TODO: is there a possibility of post-translational mutation of DNA or RNA?
 
'''

from bulbs.model import Node, Relationship
from bulbs.property import String, Integer, Float
from bulbs.utils import current_datetime


#TODO: delete all the IDs


class Meta_Protein(Node):
    element_type="Protein"
    
    Prot_ID=String(nullable=False, Index=True)
    UniprotId=String(Index=True)
    CommonNames=String(Index=True) #Names separated by ;
    Length_aa=Integer(Index=True)

class Meta_DNA(Node):
    element_type="DNA"
    
    DNA_Position=String(nullable=False, Index=True) # under the form X12341-X512341234
    Gene_name=String(Index=True)
    CommonNames=String(Index=True) #Names separated by ;
    Length_pb=Integer(Index=True)

class Meta_RNA(Node):
    element_type="RNA"
    
    RNA_ID=String(nullable=False, Index=True)
    name=String(Index=True)
    CommonNames=String(Index=True) #Names separated by ;
    Length_pb=Integer(Index=True)
    
class Meta_Complex(Node):
    element_type="Complex"
    
    Complex_ID=String(nullable=False, Index=True)
    name=String(Index=True)
    CommonNames=String(Index=True) #Names separated by ;

class Meta_SmallMolecule(Node):
    element_type="Small Molecule"

    SmallMolecule_ID=String(nullable=False,Index=True)
    ChEBI_ID=String(Index=True)
    name=String(Index=True)
    CommonNames=String(Index=True) #Names separated by;

class Fragment_CatalyticSite(Node):
    element_type="Catalytic Site"
    
    CatalyticSite_ID=String(nullable=False, Index=True)
    name=String(Index=True)
    CommonNames=String(Index=True)
    # No localization, because regroups other catalytic sites

class Fragment_Composite_CatalyticSite(Node): # Exists only within a complex, made up by several catalytic sites types
    element_type="Composit Catalytic Site"
    
    Comp_CatalyticSite_ID=String(nullable=False,Index=True)
    name=String(Index=True)
    CommonNames=String(Index=True)
    localization=String(Index=True)
    
class Fragment_ContactRegion(Node):
    element_type="Contact Region" # => Region of contact with an another protein
    
    ID=String(nullable=False, Index=True)
    name=String(Index=True)
    localization=String(Index=True)
    
class Fragment_EvoDefSite(Node): #classification introduction in a machine-readable format
    element_type="Evolutionary defined site" # => Conserved elements, evolution hotspots, etc...
    
    EvDefSite_ID=String(nullable=False, Index=True)
    name=String(Index=True)
    localization=String(Index=True)
    type=String(Index=True)

class Instantiation_PostTranslationalModificationSite(Node): # Do we need it? we can instantiate it otherwise
    element_type="Post-Translational Modification site"
    
    ID=String(nullable=False, Index=True)
    name=String(Index=True)
    Type=String(Index=True)
    localization=String(Index=True)
    
class Instantiation_Mutation(Node):
    element_type="Mutation"
    
    ID=String(nullable=False, Index=True)
    type=String(Index=True)
    Abs_Location=String(Index=True) 
    Rel_Location=String(Index=True)

class Instantiation_Epigenetic_Modification(Node): # Specific only to DNA
    element_type="Epugenetic_Modification"
    
    ID=String(nullable=False,Index=True)
    Abs_location = String(Index=True)

class Instantiation_InComplex(Node):
    element_type="In Complex"
    
    ID=String(nullable=False, Index=True)
    
class Instantiation_Localization(Node):
    element_type="Localization"
    
    ID=String(nullable=False, Index=True)
    Location=String(Index=True)
    
class Meta_Reaction(Node):
    element_type="Reaction"
    
    ID=String(nullable=False, Index=True)
    name=String(Index=True)
    type=String(Index=True)
    Free_Entalpy=String()
    Cinetic_Constant=String()

class Instance_of_Object(Node):
    element_type="Reaction"

    ID=String(nullable=False, Index=True)
    
class Instace_of_Reaction(Node):
    element_type="Instance of reaction"
    
    ID=String(nullable=False, Index=True)

class Instance_modification(Node):
    element_type="modification of instantiation"
    
    ID=String(nullable=False,Index=True)
    
class Annot_Reaction_type(Node):
    element_type="Reaction Type"
    
    ID=String(nullable=False, Index=True)
    name=String(Index=True)
    
class Annot_domain(Node):
    element_type="Domain types"
    
    name=String(nullable=False, Index=True)
    
class Annot_Locations(Node):
    element_type="Location types"
    
    name=String(nullable=False, Index=True)
    
class Annot_PostTransMod_Type(Node):
    element_type="Post Translational Modification Type"
    
    name=String(nullable=False, Index=True)

class Annot_GOTerm(Node):
    element_type="goterm"
    
    GOTermID=String(nullable=False)
    FullName=String()

class M2A_Annotates(Relationship):
    #relation between a Protein and a GOTerm
    label = "annotates"
    
    confidence=Float() # between 0 for false to 1 for sure 

#Other classes are to be constructed in relation with the Reactome data structure

class Annot_Database(Node):
    # Should not be used for the search because of a large number of nodes attached to it
    element_type="database_id"
    
    db_name=String(nullable=False, Index=True)
