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
        RelDict[ID]=copy.copy(LocList)
    Fle=file('GO.dump','w')
    pickle.dump((RelDict, SeedSet),Fle)
    return RelDict, SeedSet

def get_GO_structure(Filtr,seedList):
    '''
    Loads all of the relations between the GOs that are generalisation of the seedList GOs and that are withing the types specified in Filtr
    @param Filtr: the List of GO types we would like to get loaded into our analysis
    @type Filtr: list of strings
    
    @param seedList: the List of GO types we would like to get loaded into our analysis. It is assumed that seedList obeys the Filtr rules
    @type seedList: list of strings
    
    @return GeneralDict: ID -> Local Ontology List, Local Regulation List
    '''
    GeneralDict={}
    VisitedSet=set()
    seedList=list(seedList)
    while seedList!=[]:
        ID=seedList.pop()
        VisitedSet.add(ID)
        LocUpList=[]
        LocRegList=[]
        GONode=DatabaseGraph.GOTerm.get(ID)
        for Typ in GOUpTypes:
            generator=GONode.bothV(Typ)
            if generator!=None:
                for elt in generator:
                    subID=tr(elt).split('/')[-1][:-1]
                    if elt.Namespace in Filtr:
                        LocUpList.append(subID) 
                        if subID not in VisitedSet and subID not in seedList:
                            seedList.append(subID)
        for Typ in GORegTypes:
            generator=GONode.bothV(Typ)
            if generator!=None:
                for elt in generator:
                    subID=tr(elt).split('/')[-1][:-1]
                    if elt.Namespace in Filtr:
                        LocRegList.append(subID)
                        if subID not in VisitedSet and subID not in seedList:
                            seedList.append(subID)
        LocUpList=list(set(LocUpList))
        LocRegList=list(set(LocRegList))
        GeneralDict[ID]=(LocUpList,LocRegList)
    Fle=file('GO_structure.dump','w')
    pickle.dump(GeneralDict,Fle)
    return GeneralDict

def get_GO_Informativities():
    '''
    here calculated without any information on regulation
    '''
    init=time()
    GO_access=pickle.load(file('GO.dump','r'))
    GO_structure=pickle.load(file('GO_structure.dump','r'))
    TimesReached={}
    i=0
    l=len(GO_access)
    for key in GO_access:
        i+=1
        print 'entering',i,'/',l,time()-init
        init=time()
        toVisit=[]
        toVisit=copy.copy(GO_access[key])
        visited=[]
        while toVisit!=[]:
            elt=toVisit.pop()
            for subelt in GO_structure[elt]:
                if subelt not in visited and subelt not in toVisit:
                    toVisit.append(subelt)
            visited.append(elt)
        for elt in visited:
            if elt not in TimesReached.keys():
                TimesReached[elt]=0
            TimesReached+=1
    Fle=file('GO_Informativities.dump','w')
    pickle.dump(TimesReached,Fle)
    
            
                    
        
        