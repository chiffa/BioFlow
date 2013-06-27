'''
Created on Jun 24, 2013

@author: andrei
'''

import Graph_Declaration as GD
import copy

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

# ProtNameCollections=displayNameClusters(GD.DatabaseGraph.Protein)
# print len(ProtNameCollections)
# 
# Prot_CollNameCollections=displayNameClusters(GD.DatabaseGraph.Protein_Collection)
# print len(Prot_CollNameCollections)

def recount_dict(dico):
    correctDic={}
    for key, val in dico.items():
        if val not in correctDic.keys():
            correctDic[val]=0
        correctDic[val]+=1
    return correctDic

def recoverGenerated(generator):
    lst=[]
    if generator!=None:
        for elt in generator:
            if elt!=None:
                lst.append(str(elt).split('/')[-1][:-1])
        lst=list(set(lst))
    return lst

def updateDic(dico, lst, multiplicator):
    for elt in lst:
        if elt not in dico.keys():
            dico[elt]=0
        dico[elt]+=multiplicator

def Connected_to_Reactions(ReactionType):
    dico={}
    for Reaction in ReactionType.get_all():
        if Reaction!=None and Reaction.bothV()!=None:
            for elt in Reaction.bothV():
                ID=str(elt).split('/')[-1][:-1]
                if ID not in dico.keys():
                    dico[ID]=0 
                dico[ID]+=1
    print recount_dict(dico)
    return dico

def propagate_Connections(dico):
    dict2=copy.deepcopy(dico)
    for component_id in dico.keys():
        component_object=GD.DatabaseGraph.vertices.get(component_id)
        if component_object==None:
            print component_id
        else:
            updateDic(dict2,recoverGenerated(component_object.bothV('is_annotated')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_part_of_complex')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_localized')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_able_to_modify')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_regulant')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_part_of_collection')),dico[component_id])
    print recount_dict(dict2)
    return dict2

def merge_dictionaries(dicoList):
    result={}
    for dico in dicoList:
        for elt in dico.keys():
            if elt not in result.keys():
                result[elt]=0
            result[elt]+=dico[elt]
    return result


TR1=Connected_to_Reactions(GD.DatabaseGraph.TemplateReaction)
D1=Connected_to_Reactions(GD.DatabaseGraph.Degradation)
BR1=Connected_to_Reactions(GD.DatabaseGraph.BiochemicalReaction)


TR2=propagate_Connections(TR1)
D2=propagate_Connections(D1)
BR2=propagate_Connections(BR1)


merge1=merge_dictionaries([TR1,D1,BR1])
merge2=merge_dictionaries([TR2,D2,BR2])

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
# 
# 
# Dna_Gen=GD.DatabaseGraph.DNA.get_all()
# 
# 
# counter=0
# 
# 
# for Dna_el in Dna_Gen:
#     if Dna_el!=None and Dna_el.bothV()!=None:
#         counter+=1
#         print Dna_el.ID
#         for element in Dna_el.bothV():
#             print element
# 
# 
# 
# print counter
# 
#         
#         