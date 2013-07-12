'''
Created on Jul 11, 2013

@author: andrei
'''
from neo4j_Declarations.Graph_Declarator import DatabaseGraph
import copy
import numpy as np

edge_type_filter0=["is_same"]
edge_type_filter1_1=["is_Catalysant","is_reaction_particpant"]
edge_type_filter1_2=["is_part_of_complex","is_Regulant","is_interacting","is_Catalysant"]
edge_type_filter2=["is_part_of_collection","is_possibly_same"]
edge_type_filter3=["is_part_of_pathway","is_next_in_pathway"]

ReactionsList=[DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]

def get_Reaction_blocks():
    ReagentClusters=[]
    Seeds=set()
    for ReactionType in ReactionsList:
        for Reaction in ReactionType.get_all():
            if Reaction!=None:
                LocalList=[]
                for edge_type in edge_type_filter1_1:
                    if Reaction.bothV(edge_type)!=None:
                        for elt in Reaction.bothV(edge_type):
                            ID=str(elt).split('/')[-1][:-1]
                            LocalList.append(ID)
                            Seeds.add(ID)
                ReagentClusters.append(copy.copy(LocalList))
    return ReagentClusters, Seeds

def get_blocks_expansions(SeedSet):
    BlockClusters=[]
    for element in SeedSet:
        SeedNode=DatabaseGraph.vertices.get(element)
        LocalList=[element]
        for edge_type in edge_type_filter1_2:
            if SeedNode.bothV(edge_type)!=None:
                for elt in SeedNode.bothV(edge_type):
                    ID=str(elt).split('/')[-1][:-1]
                    LocalList.append(ID)
        BlockClusters.append(copy.copy(LocalList))
    return BlockClusters

def get_UNIPROT_equivalences(SuperSeed):  # expanded seed should be applied
    NodeEquivalences=[]
    for element in 
    
    
def get_HiNT_interactions():

def Connected_to_Reactions(ReactionType):
    dico={}
    for Reaction in ReactionType.get_all():
        if Reaction!=None:
            if Reaction.bothV("is_Catalysant")!=None:
                for elt in Reaction.bothV("is_Catalysant"):
                    ID=str(elt).split('/')[-1][:-1]
                    if ID not in dico.keys():
                        dico[ID]=0 
                    dico[ID]+=1   
            if Reaction.bothV("is_reaction_particpant")!=None:
                for elt in Reaction.bothV("is_reaction_particpant"):
                    ID=str(elt).split('/')[-1][:-1]
                    if ID not in dico.keys():
                        dico[ID]=0 
                    dico[ID]+=1
    return dico

def getMatrix():
    NodeNumber2MatrixNumber={}
    MatrixNumber2NodeNumbe={}
    ID2displayName={}
    ID2Type={}
    RelationSet=set()
    DicosList=[]
    IDSet=set()
    
    # Connect the groups of ingredients that share the same reactions
    
    # Retrieve seeds for the matrix computation
    for Reaction in ReactionsList:
        DicosList.append(Connected_to_Reactions(Reaction))
    for Dico in DicosList:
        IDSet.update(Dico.keys())
    # Expand the ID list to include all the edges likely to be reached via "
    
    # for each ID, create the forwards and reverse mappings to matrixNumbers(just increase it)
    for key in IDSet:
        
    

getMatrix()
