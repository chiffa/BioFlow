'''
Created on Jun 15, 2013

@author: andrei
'''
import logging
from bulbs.neo4jserver import Graph as Neo4jGraph
import Neo4j_typeDec_new as DDT
import DictGen as DG

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='dynamics_full.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

class Graph(Neo4jGraph):
    
    def __init__(self, config=None):
        super(Graph, self).__init__(config)
        
        #Annotations
        self.Location=self.build_proxy(DDT.Location)
        self.AnnotNode=self.build_proxy(DDT.AnnotNode)
        self.Originating_Organism=self.build_proxy(DDT.Originating_Organism)
        self.Pathway=self.build_proxy(DDT.Pathway)
        
        #Simple Compounds
        self.DNA=self.build_proxy(DDT.DNA)
        self.RNA=self.build_proxy(DDT.RNA)
        self.Protein=self.build_proxy(DDT.Protein)
        self.SmallMolecule=self.build_proxy(DDT.SmallMolecule)
        self.PhysicalEntity=self.build_proxy(DDT.PhysicalEntity)
        
        #And composite ones
        self.Complex=self.build_proxy(DDT.Complex)
        #That can be instantiated
        self.Instance=self.build_proxy(DDT.Instance)
        #By possible instantiating events
        self.ModificationFeature=self.build_proxy(DDT.ModificationFeature)
        
        #And that can be grouped as collections
        self.DNA_Collection=self.build_proxy(DDT.DNA_Collection)
        self.RNA_Collection=self.build_proxy(DDT.RNA_Collection)
        self.Protein_Collection=self.build_proxy(DDT.Protein_Collection)
        self.SmallMolecule_Collection=self.build_proxy(DDT.SmallMolecule_Collection)
        self.PhysicalEntity_Collection=self.build_proxy(DDT.PhysicalEntity_Collection)
        self.Complex_Collection=self.build_proxy(DDT.Complex_Collection)
        
        #And especially participate into reaction sets
        self.TemplateReaction=self.build_proxy(DDT.TemplateReaction)
        self.Degradation=self.build_proxy(DDT.Degradation)
        self.BiochemicalReaction=self.build_proxy(DDT.BiochemicalReaction)
        
        #Pointers from the Simple Compounds to their annotations
        self.is_localized=self.build_proxy(DDT.is_localized)
        self.is_annotated=self.build_proxy(DDT.is_annotated)
        self.is_originating_in_organism=self.build_proxy(DDT.is_originating_in_organism)
        self.is_part_of_pathway=self.build_proxy(DDT.is_part_of_pathway)
        self.is_next_in_pathway=self.build_proxy(DDT.is_next_in_pathway)
        
        #And from Complex Compounds to the simple Compounds they are made of
        self.is_part_of_complex=self.build_proxy(DDT.is_part_of_complex)
        
        #That can be instantiated
        self.is_modified_to=self.build_proxy(DDT.is_modified_to)
        #With Instantiators
        self.is_able_to_modify=self.build_proxy(DDT.is_able_to_modify)
        #Or belong to collections        
        self.is_part_of_collection=self.build_proxy(DDT.is_part_of_collection)
        
        #And contribute to reactions
        self.is_catalysant=self.build_proxy(DDT.is_catalysant)
        self.is_regulant=self.build_proxy(DDT.is_regulant) # regulates not a reaction, but a compound activity
        self.is_reaction_participant=self.build_proxy(DDT.is_reaction_participant)
        
        
#now, let's connect the graph and fill it with data from the etree parsing

DatabaseGraph=Graph()

LocalDict={} # accelerated access pointer to the objects

def InsertCellLocations():
    for Loc in DG.CellularLocations.keys():
        LocalDict[Loc]=DatabaseGraph.Location.create(ID=Loc, displayName=DG.CellularLocations[Loc])

def MinimalAnnotInsert(primary,reflist):
    for Type in reflist.keys():
        if Type!='name' and reflist[Type]!='' and reflist[Type]!=[]:
            secondary=DatabaseGraph.AnnotNode.create(ptype=Type,payload=reflist[Type])
            print primary, secondary, primary.ID
            DatabaseGraph.is_annotated.create(primary, secondary, costum_from=primary.ID,costum_to='Annotation')

def MetaInsert(function,dico):
    length=len(dico)
    counter=0
    for key in dico.keys():
        counter+=1
        print '\n',counter,'/',length,'\n'
        primary=function.create(ID=key, displayName=dico[key]['displayName'], localization=dico[key]['cellularLocation'])
        LocalDict[key]=primary
        MinimalAnnotInsert(LocalDict[key], dico[key]['references'])
        if 'cellularLocation' in dico[key].keys():
            secondary=LocalDict[dico[key]['cellularLocation']]
            DatabaseGraph.is_localized.create(primary, secondary,costum_from=primary.ID,costum_to=secondary.ID)
            # TODO add ModificationFeature insertion
        if 'modification' in dico[key].keys():
            for modification in dico[key]['modification']:
                if 'location' in modification.keys() and 'modification' in modification.keys():
                    LocMod=DatabaseGraph.ModificationFeature.create(ID = modification['ID'], type="post-translational_Mod", location=modification['location'], displayName=modification['modification'])
                    DatabaseGraph.is_able_to_modify.create(primary,LocMod,costum_from=primary.ID,costum_to=LocMod.ID)


def CollectionRefsInsert(primaryCollection):
    for key in primaryCollection.keys():
        for ref in primaryCollection[key]['collectionMembers']:
            DatabaseGraph.is_part_of_collection.create(LocalDict[key],LocalDict[ref],costum_from=LocalDict[key].ID,costum_to=LocalDict[ref].ID)

def ComplexPartsInsert():
    for key in DG.Complexes.keys():
        for part in DG.Complexes[key]['parts']:
            if 'Stoichiometry' not in part:
                DatabaseGraph.is_part_of_complex.create(LocalDict[key],LocalDict[part],costum_from=LocalDict[key].ID,costum_to=LocalDict[part].ID)

def ReactionInsert(function,dico):
    for key in dico.keys():
        LocalDict[key]=function.create(ID=key,displayName=dico[key]['displayName'])
        MinimalAnnotInsert(LocalDict[key], dico[key]['references'])
        for subkey in dico[key].keys():
            if subkey in ['left','right']:
                for elt in dico[key][subkey]:
                    DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[elt],side=subkey,costum_from=LocalDict[key].ID,costum_to=LocalDict[elt].ID) 
            if subkey=='product':
                DatabaseGraph.is_reaction_participant.create(LocalDict[key],LocalDict[dico[key][subkey]],costum_from=LocalDict[key].ID,costum_to=LocalDict[dico[key][subkey]].ID)

def CatalysisInsert():
    for key in DG.Catalysises.keys():
        if 'controller' in DG.Catalysises[key].keys() and 'controlled' in DG.Catalysises[key].keys():
            if DG.Catalysises[key]['controlled'] in LocalDict.keys() and DG.Catalysises[key]['controller'] in LocalDict.keys():
                if 'ControlType' in DG.Catalysises[key].keys():
                    primary=LocalDict[DG.Catalysises[key]['controller']]
                    secondary=LocalDict[DG.Catalysises[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary, secondary, ID=key, controlType=DG.Catalysises[key]['ControlType'],costum_from=primary.ID,costum_to=secondary.ID)
                else:
                    primary=LocalDict[DG.Catalysises[key]['controller']]
                    secondary=LocalDict[DG.Catalysises[key]['controlled']]
                    LocalDict[key]=DatabaseGraph.is_catalysant.create(primary, secondary, ID=key, controlType='UNKNOWN',costum_from=primary.ID,costum_to=secondary.ID)
            else: 
                logging.debug("\t%s : %s, %s, %s", key, DG.Catalysises[key],DG.Catalysises[key]['controlled'] in LocalDict.keys(),DG.Catalysises[key]['controller'] in LocalDict.keys())
        else: 
            logging.debug("%s : %s, %s, %s,", key, DG.Catalysises[key],'controller' in DG.Catalysises[key].keys(),'controlled' in DG.Catalysises[key].keys())



def ModulationInsert():
    for key in DG.Modulations.keys():
        primary=LocalDict[DG.Modulations[key]['controller']]
        secondary=LocalDict[DG.Modulations[key]['controlled']]
        LocalDict[key]=DatabaseGraph.is_regulant.create(primary, secondary, ID=key,controlType=DG.Modulations[key]['controlType'],costum_from=primary.ID,costum_to=secondary.ID)

def getOneMetaSet(function):
    for MetaKey in function.get_all():
        if MetaKey!=None:
            LocalDict[MetaKey.ID]=MetaKey
        
def getAllMetaSets():
    functionList=[
                  DatabaseGraph.DNA,
                  DatabaseGraph.DNA_Collection,
                  DatabaseGraph.RNA,
                  DatabaseGraph.RNA_Collection,
                  DatabaseGraph.SmallMolecule,
                  DatabaseGraph.SmallMolecule_Collection,
                  DatabaseGraph.Protein,
                  DatabaseGraph.Protein_Collection,
                  DatabaseGraph.Complex,
                  DatabaseGraph.Complex_Collection,
                  DatabaseGraph.PhysicalEntity,
                  DatabaseGraph.PhysicalEntity_Collection
                  ]
    
    for function in functionList:
        getOneMetaSet(function)

# InsertCellLocations()
# 
# MetaInsert(DatabaseGraph.DNA, DG.Dnas)
# MetaInsert(DatabaseGraph.DNA_Collection, DG.Dna_Collections)
# MetaInsert(DatabaseGraph.RNA, DG.Rnas)
# MetaInsert(DatabaseGraph.RNA_Collection, DG.Rna_Collections)
# MetaInsert(DatabaseGraph.SmallMolecule, DG.SmallMolecules)
# MetaInsert(DatabaseGraph.SmallMolecule_Collection, DG.SmallMolecule_Collections)
# MetaInsert(DatabaseGraph.Protein, DG.Proteins)
# MetaInsert(DatabaseGraph.Protein_Collection, DG.Protein_Collections)
# MetaInsert(DatabaseGraph.Complex, DG.Complexes)
# MetaInsert(DatabaseGraph.Complex_Collection, DG.Complex_Collections)
# MetaInsert(DatabaseGraph.PhysicalEntity, DG.PhysicalEntities)
# MetaInsert(DatabaseGraph.PhysicalEntity_Collection, DG.PhysicalEntity_Collections)
# CollectionRefsInsert(DG.Dna_Collections)
# CollectionRefsInsert(DG.Rna_Collections)
# CollectionRefsInsert(DG.SmallMolecule_Collections)
# CollectionRefsInsert(DG.Protein_Collections)
# CollectionRefsInsert(DG.Complex_Collections)
# CollectionRefsInsert(DG.PhysicalEntity_Collections)
#
# ComplexPartsInsert()
# 
# ## Meta insert finished
# ReactionInsert(DatabaseGraph.TemplateReaction, DG.TemplateReactions)
# ReactionInsert(DatabaseGraph.Degradation, DG.Degradations)
# ReactionInsert(DatabaseGraph.BiochemicalReaction, DG.BiochemicalReactions)
# 
# ## Reaction insert finished
# CatalysisInsert()
# ModulationInsert()