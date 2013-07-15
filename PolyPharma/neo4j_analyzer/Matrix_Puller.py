'''
Created on Jul 11, 2013

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

edge_type_filter0_1=["is_part_of_collection"]
edge_type_filter0_2=["is_same"]
edge_type_filter1_1=["is_Catalysant","is_reaction_particpant"]
edge_type_filter1_2=["is_part_of_complex", "is_Regulant"]
edge_type_filter1_3=["is_interacting"]
edge_type_filter2=["is_possibly_same"]
val1=0.9
val2=0.5
val3=0.25
DfactorDict={"Group":val1,
             "Same":val1,
             "Reaction":val2,
             "Contact_interaction":val2,
             "possibly_same":val3,
             }

# edge_type_filter3=["is_part_of_pathway","is_next_in_pathway"]

# Matrix filled with 0 (full dissipation, no connection), 1 (no dissipation, aliases) or r in [0,1] for partial dissipations



ReactionsList=[DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]

def get_Reaction_blocks():
    ReagentClusters=[]
    Seeds=set()
    count=0
    for ReactionType in ReactionsList:
        for Reaction in ReactionType.get_all():
            if Reaction!=None:
                LocalList=[]
                for edge_type in edge_type_filter1_1:
                    if Reaction.bothV(edge_type)!=None:
                        for elt in Reaction.bothV(edge_type):
                            ID=str(elt).split('/')[-1][:-1]
                            if ID not in IDFilter:
                                LocalList.append(ID)
                                Seeds.add(ID)
                                count+=1
                ReagentClusters.append(copy.copy(LocalList))
    return ReagentClusters, Seeds, count

def get_expansion(SubSeed,edge_type_filter):
    Clusters=[]
    SuperSeed=set()
    SuperSeed.update(SubSeed)
    count=0
    for element in SubSeed:
        SeedNode=DatabaseGraph.vertices.get(element)
        LocalList=[element]
        for edge_type in edge_type_filter:
            if SeedNode.bothV(edge_type)!=None:
                for elt in SeedNode.bothV(edge_type):
                    ID=str(elt).split('/')[-1][:-1]
                    if ID not in IDFilter:
                        LocalList.append(ID)
                        SuperSeed.add(ID)
                        count+=1
        if len(LocalList)>1:
            Clusters.append(copy.copy(LocalList))
    return Clusters, SuperSeed, count

def build_correspondances(IDSet):
    NodeID2MatrixNumber={}
    MatrixNumber2NodeID={}
    ID2displayName={}
    ID2Type={}
    ID2Localization={}
    counter=0
    LocationBufferDict={}
    for ID in IDSet:
        NodeID2MatrixNumber[ID]=counter
        MatrixNumber2NodeID[counter]=ID
        Vertex=DatabaseGraph.vertices.get(ID)
        ID2displayName[ID]=Vertex.displayName
        ID2Type[ID]=Vertex.element_type
        if Vertex.localization!=None:
            ID2Localization[ID]=request_location(LocationBufferDict,Vertex.localization)
        counter+=1
    return NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization

def request_location(LocationBufferDict,location):
    location=str(location)
    if location in LocationBufferDict.keys():
        return LocationBufferDict[location]
    else:
        generator=DatabaseGraph.Location.index.lookup(ID=location)
        if generator!=None:
            for elt in generator:
                LocationBufferDict[location]=str(elt.displayName)
                return str(elt.displayName)


def getMatrix(decreaseFactorDict):
    init=time()
    # Connect the groups of ingredients that share the same reactions1
    # Retrieve seeds for the matrix computation
    ReactLinks, InitSet, c = get_Reaction_blocks()
    print len(ReactLinks), len(InitSet), c
    print time()-init
    t=time()
    GroupLinks, GroupSet, c = get_expansion(InitSet,edge_type_filter0_1)
    print len(GroupLinks), len(GroupSet),c
    print time()-init, time()-t
    t=time()
    SecLinks, SecSet, c = get_expansion(GroupSet,edge_type_filter1_2)
    print len(SecLinks), len(SecSet), c
    print time()-init, time()-t
    t=time()
    UP_Links, UPSet, c = get_expansion(SecSet,edge_type_filter0_2)
    print len(UP_Links), len(UPSet), c
    print time()-init, time()-t
    t=time()
    HiNT_Links, FullSet, c = get_expansion(UPSet,edge_type_filter1_3)
    print len(HiNT_Links), len(FullSet), c
    print time()-init, time()-t
    t=time()
    Super_Links, ExpSet, c = get_expansion(FullSet,edge_type_filter2)
    print len(Super_Links), len(ExpSet), c
    print time()-init, time()-t
    t=time()
    
    loadLen=len(ExpSet)
    ValueMatrix=lil_matrix((loadLen,loadLen))
    
    # Fill in the matrix with the values
    # Take an impact vector
    # Continue multiplications as long as needed for convergence
    
    # export the matrix as a flat file
    #    => Most significantly touched elements, especially in the UNIPORT
    #    => Get the vector of affected proteins, then multiply it over the transfer
    #        Matrix until an equilibrium is reached.
    
    # Pay attention to the criticality spread => vector shoud increase exponentially for the important prots, effectively shutting down the whole system
    # But not in the case of "unimportant proteins"
    
    # => Assymetric influence matrices (causality followship)
    # Markov clustering linalgebra on sparce matrices to accelerate all this shit?
    
    # We could actually envision it as a chain reaction in a nuclear reactor, leading either to a reaction spiraling out of control (total functional shutdown, at least for a
    # given function.
    
    # Idea behind the eigenvectors: if we generate random sets of genes perturbating the network, some combination would lead to a way more powerful effect when propagated
    # in a markovian, turn-based network (runaway), whereas other sets will lead to a lighter runaway. A way to estimate runaway specifics of protein-protein interaction network
    # The strongest runaway would be generated by the highest absolute-value link
    # ACHTUNG!!!!! all the lanes have to be normalized later on to be effectively corresponding to a Markov model!!!! => Not in our case, since some proteins are affecting
    # SEVERAL other proteins at the same time
    
    # Limitations: no physical-path toxicity (such as rising pH, changing the O2 content or depleting ATP/ADP)
    
    print "building correspondances", time()-init, time()-t
    t=time()
    
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization = build_correspondances(ExpSet)
    
    
    print "building the ValMatrix", time()-init, time()-t
    t=time()

    for group in ReactLinks:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Reaction"],1-1e3)
    for group in GroupLinks:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Group"],1-1e3)
    for group in SecLinks:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1-1e3)
    for group in UP_Links:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=decreaseFactorDict["Same"]
    for group in  HiNT_Links:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1-1e3)
    for group in  Super_Links:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["possibly_same"],1-1e3)        
    
    print "entering eigenvect computation", time()-init, time()-t
    t=time()
    eigenvals, eigenvects = eigsh(ValueMatrix,1000)
    print eigenvals
    output = file('eigenvals.csv','w')
    output.write(eigenvals)
    pickleDump=file('pickleDump.dump','w')
    pickle.dump(eigenvects, pickleDump)
    pickleDump.close()
    pickleDump2=file('pickleDump2.dump','w')
    pickle.dump((NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization),pickleDump2)
    pickleDump2.close()
    pickleDump3=file('pickleDump3.dump','w')
    pickle.dump(ValueMatrix,pickleDump3)
    pickleDump3.close()
    
    print time()-init, time()-t

def get_eigenvect_Stats():
    init=time()
    pickleDump=file('pickleDump.dump','r')
    eigenvects=pickle.load(pickleDump)
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization = pickle.load(pickleDump2)
    pickleDump3=file('pickleDump3.dump','r')
    ValueMatrix=pickle.load(pickleDump3)
    eigenVectsList=[]
    SuperIndex={}
    CounterIndex={}
    print 'depicked',time()-init
    for i in range(0,eigenvects.shape[1]):
        eigenVectsList.append(eigenvects[:,i])
    i=0
    for eigenvect in eigenVectsList:
        i+=1
        normalized=np.multiply(eigenvect,eigenvect)
        for k in range(0,10):
            index=np.argmax(normalized, 0)
            if MatrixNumber2NodeID[index] not in CounterIndex.keys():
                CounterIndex[MatrixNumber2NodeID[index]]=0
                if MatrixNumber2NodeID[index] in ID2Localization.keys():
                    SuperIndex[MatrixNumber2NodeID[index]]=(ID2Type[MatrixNumber2NodeID[index]], ID2displayName[MatrixNumber2NodeID[index]], ID2Localization[MatrixNumber2NodeID[index]])
                else: 
                    SuperIndex[MatrixNumber2NodeID[index]]=(ID2Type[MatrixNumber2NodeID[index]], ID2displayName[MatrixNumber2NodeID[index]])
            CounterIndex[MatrixNumber2NodeID[index]]+=normalized[index]
            normalized[index]=0
    srtd=sorted(CounterIndex.iteritems(), key=operator.itemgetter(1),reverse=True)
    print IDFilter
    for key,val in srtd[:100]:
        print val, key, SuperIndex[key], key in IDFilter
        if key in IDFilter:
            print 'error on key: ', key
    return CounterIndex, SuperIndex
    
def processEigenVectors():
    init=time()
    pickleDump=file('pickleDump.dump','r')
    eigenvects=pickle.load(pickleDump)
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization = pickle.load(pickleDump2)
    pickleDump3=file('pickleDump3.dump','r')
    ValueMatrix=pickle.load(pickleDump3)
    eigenVectsList=[]
    SuperIndex=[]
    print 'depicked',time()-init
    for i in range(0,eigenvects.shape[1]):
        eigenVectsList.append(eigenvects[:,i])
    i=0
    print len(eigenVectsList)
    for eigenvect in eigenVectsList:
        i+=1
        normalized=np.multiply(eigenvect,eigenvect)
        indexes={}
        for k in range(0,10):
            index=np.argmax(normalized, 0)
            print eigenvect.shape
            print normalized.shape
            indexes[MatrixNumber2NodeID[index]]=(normalized[index], ID2Type[MatrixNumber2NodeID[index]], ID2displayName[MatrixNumber2NodeID[index]])
            normalized[index]=0
        SuperIndex.append(indexes)
    for subIndex in SuperIndex:
        for key in subIndex.keys():
            print '\t', key, subIndex[key]
        print '<=========================>'
    return SuperIndex

getMatrix(DfactorDict)

# processEigenVectors()

get_eigenvect_Stats()

# TODO: create GO and Pathway Structure access
# Calibrate the values so that after ~3 transitions the correlation vanishes on average (Follow Pamela Silver Approach)
