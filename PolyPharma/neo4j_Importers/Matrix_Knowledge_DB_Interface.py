__author__ = 'ank'
"""
Contains all the tools necessary to map GO ontology and Pathway classification from the database to an Adjacency and
Laplacian graph.
"""

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
import copy
from time import time
import pickle
import operator
from random import shuffle
import math
from scipy.stats.kde import  gaussian_kde
import numpy as np
from pylab import plot, hist, show
from PolyPharma.neo4j_analyzer.Matrix_Interactome_DB_interface import MatrixGetter
from PolyPharma.configs import Dumps
from PolyPharma.neo4j_analyzer.knowledge_access import acceleratedInsert

# Creates an instance of MatrixGetter and loads pre-computed values
MG = MatrixGetter(True, False)
MG.fast_load()

# specify the relations that lead to a more general or to an equally regulated node.
GOUpTypes = ["is_a_go", "is_part_of_go"]
GORegTypes = ["is_Regulant"]


class GO_Interface(object):
    """
    General calss to recover all the informations associated with GO from database and buffer them for further use.
    """

    def __init__(self):
        pass

    def get_GO_access(self, Filtr):
        '''
        Loads all of the relations between the UNIPROTs and GOs as one giant dictionary

        :param Filtr: the List of GO types we would like to get loaded into our analysis
        :type Filtr: list of strings

        :return: RelDict, NodeID -> list of IDs of GO annotations nodes associated to it
        :return: SeedSet, List of GO nodes reached from the NodeID
        '''
        RelDict = {}
        SeedSet = set()
        UPs_without_GO = 0
        for ID in MG.Uniprots:
            UP_Specific_GOs = []
            Root = DatabaseGraph.UNIPORT.get(ID)
            Node_gen = Root.bothV("is_go_annotation")
            if Node_gen:
                for GO in Node_gen:
                    if GO.Namespace in Filtr:
                        GOID = str(GO).split('/')[-1][:-1]
                        UP_Specific_GOs.append(GOID)
                        SeedSet.add(GOID)
            if UP_Specific_GOs == []:
                UPs_without_GO += 1
                print "Warning: UP without GO has been found. Database ID: %s, name: %s!!!!!!" % (ID, MG.ID2displayName[ID])
            else:
                RelDict[ID] = copy.copy(UP_Specific_GOs)
        print UPs_without_GO
        pickle.dump((RelDict, SeedSet), file(Dumps.GO_dump, 'w'))  #TODO": correct dumping here
        return RelDict, SeedSet

    def get_GO_structure(self, Filtr,seedSet):
        '''
        Loads all of the relations between the GOs that are generalisation of the seedList GOs and that are withing the types specified in Filtr

        :param Filtr: the List of GO types we would like to get loaded into our analysis
        :type Filtr: list of strings
        :param seedSet: set of GO types we would like to get loaded into our analysis. It is assumed that seedList obeys the Filtr rules
        :type seedSet: set of strings

        :return: GeneralDict, ID -> Local Ontology List, Local Regulation List
        :rtype: dict
        '''
        GeneralDict = {}
        VisitedSet = set()
        seedList = list(seedSet)
        GO_Names = {}
        GO_IDs = {}
        rev_GO_IDs = {}
        while seedList:
            ID = seedList.pop()
            VisitedSet.add(ID)
            LocUpList = []
            LocRegList = []
            GONode = DatabaseGraph.GOTerm.get(ID)
            GO_Names[ID] = str(GONode.displayName)
            GO_IDs[ID] = str(GONode.ID)
            rev_GO_IDs[str(GONode.ID)] = ID
            for Typ in GOUpTypes:
                generator = GONode.outV(Typ)
                if generator:
                    for elt in generator:
                        subID = str(elt).split('/')[-1][:-1]
                        if elt.Namespace in Filtr:
                            LocUpList.append(subID)
                            if subID not in VisitedSet and subID not in seedList:
                                seedList.append(subID)
            for Typ in GORegTypes:
                generator = GONode.outV(Typ)
                if generator:
                    for elt in generator:
                        subID = str(elt).split('/')[-1][:-1]
                        if elt.Namespace in Filtr:
                            LocRegList.append(subID)
                            if subID not in VisitedSet and subID not in seedList:
                                seedList.append(subID)
            LocUpList = list(set(LocUpList))
            LocRegList = list(set(LocRegList))
            GeneralDict[ID] = (LocUpList, LocRegList)
        Fle = file('GO_structure.dump','w')
        pickle.dump(GeneralDict,Fle)  #TODO": correct dumping here
        Fle2 = file('GO_names.dump','w')
        pickle.dump(GO_Names,Fle2)  #TODO": correct dumping here
        Fle3 = file('GO_IDs.dump','w')
        pickle.dump(GO_IDs,Fle3)  #TODO": correct dumping here
        Fle4 = file('rev_GO_IDs.dump','w')
        pickle.dump(rev_GO_IDs,Fle4)  #TODO": correct dumping here
        return GeneralDict


    def get_GO_Informativities(self):
        '''
        here calculated without any information on regulation
        '''
        init=time()
        GO_access=pickle.load(file('GO.dump','r'))[0]  #TODO": correct dumping here
        GO_structure=pickle.load(file('GO_structure.dump','r'))  #TODO": correct dumping here
        TimesReached={}
        i=0
        l=len(GO_access)
        accelerationDict={}
        Reverse_Dict={}
        for key in GO_access.keys():
            i+=1
            print 'entering',float(i)/float(l),time()-init
            init=time()
            toVisit=copy.copy(GO_access[key])
            visited=[]
            while toVisit:
                elt=toVisit.pop()
                InStack=[]
                vs=acceleratedInsert(GO_structure, accelerationDict, elt)
                visited=visited+vs
            visited=list(set(visited))
            for elt in visited:
                if elt not in TimesReached.keys():
                    Reverse_Dict[elt]=[]
                    TimesReached[elt]=0
                Reverse_Dict[elt].append(key)
                TimesReached[elt]+=1
        Fle=file('GO_Informativities.dump','w')
        pickle.dump(TimesReached,Fle)  #TODO": correct dumping here
        Fle2=file('accDict.dump','w')
        pickle.dump(accelerationDict,Fle2)  #TODO": correct dumping here
        Fle3=file('Reverse_dict.dump','w')
        pickle.dump(Reverse_Dict, Fle3)  #TODO": correct dumping here


    def compute_UniprotDict(self):
        """


        :return:
        """
        UniprotDict = {}
        for elt in MG.Uniprots:
            node = DatabaseGraph.UNIPORT.get(elt)
            altID = node.ID
            UniprotDict[altID] = (elt, MG.ID2displayName[elt])
            UniprotDict[elt] = altID
        pickle.dump(UniprotDict, file(Dumps.Up_dict_dump,'w'))
        return UniprotDict



    def rebuild(self):
        filtr=['biological_process']
        RelDict, SeedSet = self.get_GO_access(filtr)
        self.get_GO_structure(filtr, SeedSet)
        self.get_GO_Informativities()
        self.compute_UniprotDict()