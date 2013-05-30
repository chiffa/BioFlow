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
 TODO: is there a possibility of post-translational modification of DNA or RNA?
 
 All the information regarding the 
'''

from bulbs.model import Node, Relationship
from bulbs.property import String, Integer, Float, Dictionary, List
from bulbs.utils import current_datetime


#CancelledTODO: delete all the IDs => We need to search somehow with external IDs
#TODO: Transform database pointers into the veretexes pointing towards the nodes to allow bulbs

# TODO: pointers towards common names, common names being indexed

# TODO: add a SQL database allowing the mapping of Accession numbers and Co to the Objects?
# WE ARE NOT INTO OPTIMIZATION YET!!!

# TODO: add fragement to fragment composition relation

class Meta(Node):
    ID=String(nullable=False)
    name=String(Index=True)
    custom=String(Index=True)
    Size=String()
    
class Fragment(Node):
    ApproximateStart=Integer(nullable=False, Index=True)
    ApproxinateEnd=Integer(nullable=False, Index=True)
    location=String(nullable=False,) #Location is relative for everything except DNA fragmements
    name=String(Index=True)
    custom=String(Index=True)

class Instantiation(Node):
    ID=String(nullable=False,Index=True)
    name=String(Index=True)
    custom=String(Index=True)
    
class Reaction(Node):
    ID=String(nullable=False, Index=True)
    name=String(Index=True)
    type=String(Index=True)
    Free_Entalpy=String()
    Cinetic_Constant=String()
    frequency=Float()
    custom=String(Index=True)

#<=======================================================================================================>
#<=======================================================================================================>
#<=======================================================================================================>

class Meta_Protein(Meta):
    element_type="Protein"
    
    UniprotId=String(Index=True)
    CommonNames=List(Index=True) # No, this doesn't allow efficient indexing

class Meta_Chromosome(Meta):
    element_type="Chromosome or Chromosome fragment"

#In theory DNA is always a fragment, not an entity. Chromosome or chromosome fragment / recombination is

class Meta_RNA(Meta):
    element_type="RNA"
    
    EMBL_RNA_ID=String(Index=True)
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing
    
class Meta_Complex(Meta):
    element_type="Complex"
    
    PDB_ID=String(Index=True)
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing

class Meta_SmallMolecule(Meta):
    element_type="Small Molecule"

    ChEBI_ID=String(Index=True)
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing

#<=======================================================================================================>

class Fragment_DNA(Fragment):
    element_type="DNA Fragment"
    
    Gene_name=String(Index=True)
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing
    Length_pb=Integer(Index=True)

class Fragment_CatalyticSite(Fragment):
    element_type="Catalytic Site"
    
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing
    # No localization, because regroups other catalytic sites

class Fragment_Composit_CatalyticSite(Fragment): # Exists only within a complex, made up by several catalytic sites types
    element_type="Composit Catalytic Site"
    
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing
    
class Fragment_ContactRegion(Fragment):
    element_type="Contact Region" # => Region of contact with an another protein
    
    CommonNames=String(Index=True) # No, this doesn't allow efficient indexing
    
class Fragment_EvoDefSite(Fragment): #classification introduction in a machine-readable format
    element_type="Evolutionary defined site" # => Conserved elements, evolution hotspots, etc...
    
    EvoType=String(Index=True)

class Fragment_PTM_Site(Fragment):
    element_type="Post-Translational Modification site"

class Fragment_Mutation_Site(Fragment):
    element_type="Mutation Site"

class Fragment_Domain(Fragment):
    element_type="Domain"

    DomainType=String(Index=True)
    
#<======================================================================================================>

# Instantiation always occurs at a place defined by a Fragment

class Instantiation_PostTranslationalModification(Instantiation): # Do we need it? => Yes, for annotation purposes
    element_type="Post-Translational Modification"
    
    Type=String(Index=True)
    
class Instantiation_Mutation(Node):
    element_type="Mutation"
    
    type=String(Index=True)

class Instantiation_Epigenetic_Modification(Node): # Acts on DNA fragments
    element_type="Epugenetic_Modification"

class Instantiation_InComplex(Node):
    element_type="In Complex"

class Instantiation_Conformation(Node):
    element_type="Conformation"
    
class Instantiation_Localization(Node):
    element_type="Localization"

    Location=String(Index=True)

#<========================================================================================================>
    
class Instance_of_Meta(Meta): # Points towards original Meta and the Instantiators
    element_type="Instance of Meta"
    
    UID=String()
    
class Instance_of_Fragment(Meta): # Points towards original Meta and the Instantiators
    element_type="Instance of Meta"
    
    UID=String()
    
class Instace_of_Reaction(Reaction): # Points towards original Reaction and the Instantiators
    element_type="Instance of reaction"

    UID=String()
#<========================================================================================================>

class Instance_modification(Reaction):
    element_type="Instance change Reaction"  # Can change only the instantiation of Meta. Reaction will flow anywhere the ingredients required for it are present

#<========================================================================================================>

class Annot_Reaction_type(Node):
    element_type="Reaction type"
    
    name=String(Index=True)
    
class Annot_domain_type(Node):
    element_type="Domain type"
    
    name=String(nullable=False, Index=True)
    
class Annot_Location_type(Node):
    element_type="Location type"
    
    name=String(nullable=False, Index=True)
    
class Annot_PostTransMod_type(Node):
    element_type="Post Translational Modification Type"
    
    name=String(nullable=False, Index=True)

class Annot_GOTerm(Node):
    element_type="GO term"
    
    GOTermID=String(nullable=False)
    FullName=String()

class Annot_Database(Node):
    # Should not be used for the search because of a large number of nodes attached to it
    element_type="database name"
    
    db_name=String(nullable=False, Index=True)

class Annot_CommonName(Node):
    element_type="common name"
    
    common_name=String(nullable=False, Index=True)
    
class Annot_Evidence(Node):
    # Annotation evidence
    element_type="evidence types"
    
    evidence_type=String()
    
class Annot_EvidenceInstance(Node): # Publication or a part of publication
    element_type="evidence instance"
    
    type=String(String=True)
    evidence_instance=Dictionary()
    
class Annot_Pathway(Node):
    # A pathway
    element_type="Pathway"

#<========================================================================================================>
#<========================================================================================================>
#<========================================================================================================>

class CostumRelationship(Relationship):
    confidence=Float() # between 0 for false to 1 for sure 
    costum=String()

class R2R_NextStepInPathway(CostumRelationship):
    label="next step"

class R_A2A_InPatway(CostumRelationship):
    label="In pathway"

class toA_CommonName(Relationship):
    label="Commonly named"

class M_F2A_Annotates(CostumRelationship):
    # Relation between a meta-object and a GOTerm
    label = "annotates"

class M_F_I2A_Xref(Relationship):
    # Relationship between a meta-objects and their representations in other databases
    label = "Xref"
    
    id=String(nullable=False,Index=True)

class M2R_Participates_in_reaction(CostumRelationship):
    # Relation between a meta-object and a reaction:
    label="Participates"
    
    side=String() # left or right
    stocheometry=Integer()

class M2M_Encodes(Relationship):
    #relation between a protein and a rna/dna fragment that encodes it
    label="Encodes"

class M2R_Regulates(CostumRelationship):
    label="regulates"

class M2F_Maps_To(CostumRelationship):
    # Maps a meta-objetct to the fragment of DNA/RNA that encoded it
    label="Maps"

class M_F2M_Part_of(Relationship):
    #Relations of appartenance between Fragments and Meta-objects and Complex and other Meta-objects
    label="Part_of"
    
class M_F2IM_IF_Instance_of(CostumRelationship):
    label="Instance of"
    
class ToA_Evidence(CostumRelationship):
    label="Evidence"
    
class ToA_Evidence_type(Relationship):
    label="Evidence type"

class ToA_Typing(CostumRelationship):
    label="Types"

class IM_IF_IR2I_Instantiated_by(CostumRelationship):
    label="Instantiated by"
    
class F2M_BelongsTo(CostumRelationship):
    label="Belongs To"