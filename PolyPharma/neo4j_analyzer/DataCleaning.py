'''
Created on Jun 24, 2013

@author: andrei
'''

from neo4j_analyzer import Graph_Declaration as GD

def displayNameClusters(func):
    displayNameClusters={}
    
    for Meta in func.get_all():
        if Meta.displayName not in displayNameClusters.keys():
            displayNameClusters[Meta.displayName]=[]
        displayNameClusters[Meta.displayName].append(Meta)
    
    for clusterKey in displayNameClusters.keys():
        if len(displayNameClusters[clusterKey])<2:
            del displayNameClusters[clusterKey]
    
    return displayNameClusters

ProtNameCollections=displayNameClusters(GD.DatabaseGraph.Protein)
print len(ProtNameCollections)

Prot_CollNameCollections=displayNameClusters(GD.DatabaseGraph.Protein_Collection)
print len(Prot_CollNameCollections)

def Connected_to_Reactions(ReactionType):
    lst=[]
    for Reaction in ReactionType.get_all():
        lst.append(Reaction.outV('is_reaction_participant'))
    
    for Nde in lst:
        # if this is a complex
        
        # if this is an annotation
        
        # if this is a 
        
    
        
        