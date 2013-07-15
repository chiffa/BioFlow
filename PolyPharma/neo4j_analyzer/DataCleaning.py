'''
Created on Jun 24, 2013

@author: andrei
'''

from neo4j_Declarations.Graph_Declarator import DatabaseGraph
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


def crossLinkClusters(dictionary):
    for commonDispName in dictionary.keys():
        print 'enterCopy'
        lst=copy.deepcopy(dictionary[commonDispName])
        while lst!=[]:
            elt=lst.pop()
            print 'popped ', elt
            for elt2 in lst:
                print '\t interated over', elt2
                DatabaseGraph.is_possiby_same.create(elt,elt2)
        


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

def resume_dict(dico):
    resumeDic={}
    for key, val in dico.items():
        if val not in resumeDic.keys():
            resumeDic[val]=[]
        resumeDic[val].append(key)
    return resumeDic

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

def updateDic_experimental(dico, lst, multiplicator):
    for elt in lst:
        if elt in dico.keys():
            dico[elt]+=multiplicator

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
    print recount_dict(dico)
    return dico


# TODO: perform a join over the displayNames and eliminate elements of Collections that are not encountered in the other reactions
def propagate_Connections(dico):
    dict2=copy.deepcopy(dico)
    for component_id in dico.keys():
        component_object=DatabaseGraph.vertices.get(component_id)
        if component_object==None:
            print component_id
        else:
            #updateDic(dict2,recoverGenerated(component_object.bothV('is_annotated')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_part_of_complex')),dico[component_id])
            #updateDic(dict2,recoverGenerated(component_object.bothV('is_localized')),dico[component_id])
            #updateDic(dict2,recoverGenerated(component_object.bothV('is_able_to_modify')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_Regulant')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_part_of_collection')),dico[component_id])
            updateDic(dict2,recoverGenerated(component_object.bothV('is_possibly_same')),dico[component_id])
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



TR1=Connected_to_Reactions(DatabaseGraph.TemplateReaction)
D1=Connected_to_Reactions(DatabaseGraph.Degradation)
BR1=Connected_to_Reactions(DatabaseGraph.BiochemicalReaction)
 
 
merge1=merge_dictionaries([TR1,D1,BR1])
merge2=propagate_Connections(merge1)

# TODO: pull the list of elements attended from the most points
