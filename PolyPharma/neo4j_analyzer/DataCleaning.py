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
        if Reaction!=None and Reaction.bothV()!=None:
            for elt in Reaction.bothV():
                lst.append(str(elt).split('/')[-1][:-1])
    return len(lst), lst


# def Connected_to_Reactions(ReactionType):
#     lst=[]
#     for Reaction in ReactionType.get_all():
#         lst.append(Reaction.outV('is_reaction_participant'))
#     lst=list(set(lst))
#     lst2=[]
#     for Nde in lst:
#         if Nde!=None:
#             print Nde
#             # if this is a complex, get all it's parts
#             for elt in Nde.bothV('is_part_of_complex'):
#                 lst2.append(elt)
#             # if this is an annotation, keep it
#             for elt in Nde.outV('is_localized'):
#                 lst2.append(elt)
#             for elt in Nde.inV('is_localized'):
#                 lst2.append(elt)
#             # if this is a modification, keep it
#             for elt in Nde.outV('is_able_to_modify'):
#                 lst2.append(elt)
#             for elt in Nde.inV('is_able_to_modify'):
#                 lst2.append(elt)
#             # if it is a regulant or regulated, get it's part
#             for elt in Nde.outV('is_regulant'):
#                 lst2.append(elt)
#             for elt in Nde.inV('is_regulant'):
#                 lst2.append(elt) 
#             # if it's a collection, get all of it's elements
#             for elt in Nde.outV('is_part_of_collection'):
#                 lst2.append(elt)
#             for elt in Nde.inV('is_part_of_collection'):
#                 lst2.append(elt) 
#     lst2=list(set(lst))
#     return lst, lst2

# TODO: create a pullMatrix class


# Recover the Ids of the nodes in a form that is easily accessible to the neo4j database
# reclist.append(str(ModObj).split('/')[-1][:-1])


Dna_Gen=GD.DatabaseGraph.DNA.get_all()


counter=0


for Dna_el in Dna_Gen:
    if Dna_el!=None and Dna_el.bothV()!=None:
        counter+=1
        print Dna_el.ID
        for element in Dna_el.bothV():
            print element



print counter

        
        