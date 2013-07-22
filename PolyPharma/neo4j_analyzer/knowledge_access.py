'''
Created on Jul 16, 2013

@author: andrei
'''

from neo4j_Declarations.Graph_Declarator import DatabaseGraph
from configs import IDFilter
import copy
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import itertools
from time import time
import pickle
import numpy as np
import operator

GOUpTypes=["is_a_go","is_part_of_go"]
GORegTypes=["is_Regulant"]


def import_TargetMappings():
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    pickleDump2.close()
    return NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots

def import_RelMatrix():
    pickleDump3 = file('pickleDump3.dump','r')
    ValueMatrix = pickle.load(pickleDump3)
    pickleDump3.close()
    return ValueMatrix

def get_GO_access(Filtr):
    '''
    Loads all of the relations between the UNIPROTs and GOs as one giant dictionary
    @param Filtr: the List of GO types we would like to get loaded into our analysis
    @type Filtr: list of strings
    
    @return RelDict: NodeID -> list of IDs of GO annotations nodes associated to it
    @return SeedSet: List of GO nodes reached from the NodeID
    '''
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = import_TargetMappings()
    RelDict={}
    SeedSet=set()
    i=0
    for ID in Uniprots:
        LocList=[]
        Root=DatabaseGraph.UNIPORT.get(ID)
        generator = Root.bothV("is_go_annotation")
        if generator!=None:
            for GO in generator:
                if GO.Namespace in Filtr:
                    GOID=str(GO).split('/')[-1][:-1]
                    LocList.append(GOID)
                    SeedSet.add(GOID)
        if LocList==[]:
            i+=1
            print "error!!!!!!", ID, ID2displayName[ID]
        RelDict[ID]=copy.copy(LocList)
    print i
    Fle=file('GO.dump','w')
    pickle.dump((RelDict, SeedSet),Fle)
    return RelDict, SeedSet

def get_GO_structure(Filtr,seedSet):
    '''
    Loads all of the relations between the GOs that are generalisation of the seedList GOs and that are withing the types specified in Filtr
    @param Filtr: the List of GO types we would like to get loaded into our analysis
    @type Filtr: list of strings
    
    @param seedSet: set of GO types we would like to get loaded into our analysis. It is assumed that seedList obeys the Filtr rules
    @type seedSet: set of strings
    
    @return GeneralDict: ID -> Local Ontology List, Local Regulation List
    '''
    GeneralDict={}
    VisitedSet=set()
    seedList=list(seedSet)
    GO_Names={}
    while seedList!=[]:
        ID=seedList.pop()
        VisitedSet.add(ID)
        LocUpList=[]
        LocRegList=[]
        GONode=DatabaseGraph.GOTerm.get(ID)
        GO_Names[ID]=str(GONode.displayName)
        for Typ in GOUpTypes:
            generator=GONode.outV(Typ)
            if generator!=None:
                for elt in generator:
                    subID=str(elt).split('/')[-1][:-1]
                    if elt.Namespace in Filtr:
                        LocUpList.append(subID) 
                        if subID not in VisitedSet and subID not in seedList:
                            seedList.append(subID)
        for Typ in GORegTypes:
            generator=GONode.outV(Typ)
            if generator!=None:
                for elt in generator:
                    subID=str(elt).split('/')[-1][:-1]
                    if elt.Namespace in Filtr:
                        LocRegList.append(subID)
                        if subID not in VisitedSet and subID not in seedList:
                            seedList.append(subID)
        LocUpList=list(set(LocUpList))
        LocRegList=list(set(LocRegList))
        GeneralDict[ID]=(LocUpList,LocRegList)
    Fle=file('GO_structure.dump','w')
    pickle.dump(GeneralDict,Fle)
    Fle2=file('GO_names.dump','w')
    pickle.dump(GO_Names,Fle2)
    return GeneralDict

def get_GO_Informativities():
    '''
    here calculated without any information on regulation
    '''
    init=time()
    GO_access=pickle.load(file('GO.dump','r'))[0]
    GO_structure=pickle.load(file('GO_structure.dump','r'))
    TimesReached={}
    i=0
    l=len(GO_access)
    accelerationDict={}
    for key in GO_access.keys():
        i+=1
        print 'entering',float(i)/float(l),time()-init
        init=time()
        toVisit=[]
        toVisit=copy.copy(GO_access[key])
        visited=[]
        while toVisit!=[]:
            elt=toVisit.pop()
            InStack=[]
            vs=acceleratedInsert(GO_structure, accelerationDict, elt, InStack)
            visited.append(elt)
            visited=visited+vs
            visited=list(set(visited))
        for elt in visited:
            if elt not in TimesReached.keys():
                TimesReached[elt]=0
            TimesReached[elt]+=1
    Fle=file('GO_Informativities.dump','w')
    pickle.dump(TimesReached,Fle)

def load_GO_informativities():
    Fle = file('GO_Informativities.dump','r')
    GO_Infos = pickle.load(Fle)
    srted = sorted(GO_Infos.iteritems(), key=operator.itemgetter(1), reverse=True)
    GO_2_Names=pickle.load(file('GO_names.dump','r'))
    i=0
    for key, val in srted[:500]:
        i+=1
        print i, key, val, GO_2_Names[key]
    return GO_Infos

def load_GO_Structre():
    GO_structure = pickle.load(file('GO_structure.dump','r'))
    return GO_structure

def acceleratedInsert(GO_structure,accelerationDict, name, InStack):
    GO_2_Names=pickle.load(file('GO_names.dump','r'))
    if name in accelerationDict.keys():
        return accelerationDict[name]
    else:
        cumset=set()
        for subelt in GO_structure[name][0]:
            if subelt not in InStack:
                InStack.append(subelt)
                cumset.update(acceleratedInsert(GO_structure,accelerationDict, subelt, InStack))
        accelerationDict[name]=list(cumset)
        return accelerationDict[name]

def unzip(Supset, GO_2_Names):
    nameList=[]
    for elt in Supset:
        nameList.append(GO_2_Names[elt])
    return nameList
        
    
init=time()
# filtr=['biological_process']
# RelDict, SeedSet=get_GO_access(filtr)
# print 1, time()-init
# get_GO_structure(filtr,SeedSet)
# print 2, time()-init
# get_GO_Informativities()
load_GO_Structure()
load_GO_informativities()
print 3, time()-init

# TODO: check that each unipront points to a GO
# => Only about 100 uniprots out of 4000 do not point towards the 