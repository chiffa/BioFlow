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
from os import listdir
import random

edge_type_filter0_1=["is_part_of_collection"]
edge_type_filter0_2=["is_same"]
edge_type_filter1_1=["is_Catalysant","is_reaction_particpant"]
edge_type_filter1_2=["is_part_of_complex", "is_Regulant"]
edge_type_filter1_3=["is_interacting"]
edge_type_filter2=["is_possibly_same"]
val1=0.5
val2=0.33
val3=0.10
DfactorDict={"Group":val1,
             "Same":1,
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
    Clusters={}
    SuperSeed=set()
    SuperSeed.update(SubSeed)
    count=0
    for element in SubSeed:
        SeedNode=DatabaseGraph.vertices.get(element)
        LocalList=[]
        for edge_type in edge_type_filter:
            if SeedNode.bothV(edge_type)!=None:
                for elt in SeedNode.bothV(edge_type):
                    ID=str(elt).split('/')[-1][:-1]
                    if ID not in IDFilter:
                        LocalList.append(ID)
                        SuperSeed.add(ID)
                        count+=1
        if len(LocalList)>1:
            Clusters[element]=copy.copy(LocalList)
    return Clusters, SuperSeed, count

def build_correspondances(IDSet):
    NodeID2MatrixNumber={}
    MatrixNumber2NodeID={}
    ID2displayName={}
    ID2Type={}
    Uniprots=[]
    ID2Localization={}
    counter=0
    LocationBufferDict={}
    for ID in IDSet:
        NodeID2MatrixNumber[ID]=counter
        MatrixNumber2NodeID[counter]=ID
        Vertex=DatabaseGraph.vertices.get(ID)
        ID2displayName[ID]=Vertex.displayName
        ID2Type[ID]=Vertex.element_type
        if Vertex.element_type=="UNIPROT":
            Uniprots.append(ID)
        if Vertex.localization!=None:
            ID2Localization[ID]=request_location(LocationBufferDict,Vertex.localization)
        counter+=1
    return NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots

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


def getMatrix(decreaseFactorDict, numberEigvals):
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
    
    # If a specific set of GO_Terms is put down, we can say that the function they describe is down.
    # Recall v.s. precision for a GO array for a perturbed protein set?
    # Non-randomness of a recall?
    # Pathway structure?
    
    # Method extendable to inhibition / activation binaries, by introducing positive / negative values for the matrix
    
    print "building correspondances", time()-init, time()-t
    t=time()
    
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = build_correspondances(ExpSet)
    
    
    print "building the ValMatrix", time()-init, time()-t
    t=time()
    
    # Group node definintion have to be corrected so they are not all related together but instead are linked towards the central "group" node!!!!
    
    # TODO: update all this element so that it incorporates the new linkage 
    for group in ReactLinks:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Reaction"],1)
    for key in GroupLinks.keys():
        for val in GroupLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Group"],1)
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Group"],1)
    for key in SecLinks.keys():
        for val in SecLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
    for key in UP_Links.keys():
        for val in UP_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Same"],1)
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Same"],1)
    for key in HiNT_Links.keys():
        for val in HiNT_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
    for key in Super_Links.keys():
        for val in Super_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["possibly_same"],1)
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["possibly_same"],1)    

    print "entering eigenvect computation", time()-init, time()-t
    t=time()
    eigenvals, eigenvects = eigsh(ValueMatrix,numberEigvals)
    print eigenvals
    output = file('eigenvals.csv','w')
    output.write(eigenvals)
    pickleDump=file('pickleDump.dump','w')
    pickle.dump((eigenvals, eigenvects), pickleDump)
    pickleDump.close()
    pickleDump2=file('pickleDump2.dump','w')
    pickle.dump((NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots),pickleDump2)
    pickleDump2.close()
    pickleDump3=file('pickleDump3.dump','w')
    pickle.dump(ValueMatrix,pickleDump3)
    pickleDump3.close()
    
    print time()-init, time()-t

def get_Descriptor(MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, index):
    if MatrixNumber2NodeID[index] in ID2Localization.keys():
        return (ID2Type[MatrixNumber2NodeID[index]], ID2displayName[MatrixNumber2NodeID[index]], ID2Localization[MatrixNumber2NodeID[index]])
    else: 
        return (ID2Type[MatrixNumber2NodeID[index]], ID2displayName[MatrixNumber2NodeID[index]])
    
def get_eigenvect_Stats():
    init=time()
    pickleDump=file('pickleDump.dump','r')
    eigenvals, eigenvects = pickle.load(pickleDump)
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
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
                SuperIndex[MatrixNumber2NodeID[index]]=get_Descriptor(MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, index)
            CounterIndex[MatrixNumber2NodeID[index]]+=normalized[index]
            normalized[index]=0
    srtd=sorted(CounterIndex.iteritems(), key=operator.itemgetter(1),reverse=True)
    print IDFilter
    for key,val in srtd[:100]:
        print val, key, SuperIndex[key], key in IDFilter
        if key in IDFilter:
            print 'error on key: ', key
    print eigenvals
    return CounterIndex, SuperIndex

def UniprotCalibrate(rounds,depth, filename, Rdom):
    '''
    Checks if the decrease corresponds on average to the value predicted for natural networks by 
    P. Silver in E.Coli. One single propagation iteration
    # NOTICE: it might be better to perform the information collection only for the other uniprot- proteins
    '''
    init=time()
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    pickleDump3=file('pickleDump3.dump','r')
    ValueMatrix=pickle.load(pickleDump3)
    print 'unpickled in:', time()-init
    Finale=MainUprotCalLoop(Uniprots,ValueMatrix,NodeID2MatrixNumber,rounds,MatrixNumber2NodeID,ID2displayName,ID2Type,ID2Localization,depth,Rdom)
    write=file(filename, 'w')
    pickle.dump(Finale, write)
    write.close()
    print 'round completed in:', time()-init      


def MainUprotCalLoop(Uniprots,ValueMatrix,NodeID2MatrixNumber,rounds,MatrixNumber2NodeID,ID2displayName,ID2Type,ID2Localization,depth,Rdom):
    iterations=1
    ReUniprots=copy.copy(Uniprots)
    Finale=[]
    if Rdom:
        iterations=3
    for i in range(0,iterations):
        Subfinale=[]
        if Rdom:
            Uniprots=random.sample(Uniprots, 200)
        for ID in Uniprots:
            Vector=np.zeros((ValueMatrix.shape[1],1))
            Vector[NodeID2MatrixNumber[ID]]=1.0
            ForbidList=[]
            for i in range(0,rounds):
                for k in Vector.nonzero():
                    ForbidList.append(k)
                Vector=ValueMatrix*Vector
                for k in ForbidList:
                    Vector[k]=0.0
            LocalSample={ID:get_Descriptor(MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, NodeID2MatrixNumber[ID])}
            nzIndexes = Vector.nonzero()[0]
            redepth=min(depth,len(nzIndexes))
            suplist=random.sample(nzIndexes,redepth)
            for index in suplist:
                index=int(index)
                value=0
                if np.linalg.norm(Vector,1) < 10e-3:
                    print 'error', ID, np.linalg.norm(Vector,1)
                    value='infinity'
                else:
                    value=Vector[index]/np.linalg.norm(Vector,1)
                LocalSample[MatrixNumber2NodeID[index]]=(Vector[index], value, get_Descriptor(MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, index))
            Subfinale.append(LocalSample)
        Finale.append(Subfinale)
    return Finale

def mass_Calibrate(maxrange,depth,Rdom=False):
    '''
    Performs several rounds of calibration, actually recomputing the figure presented by P. Silver
    '''
    for i in range(2, maxrange+1):
        filename=str('calibrate'+str(i)+'.dump')
        UniprotCalibrate(i,depth,filename,Rdom)
    
def treat_Calibration():
    '''
    Performs the analysis of calibration
    '''
    DicList={}
    FnameList=[]
    Filenames = listdir('.')
    for filename in Filenames:
        if 'calibrate' in filename:
            FnameList.append(filename)
    for fname in FnameList:
        DicList[fname.split('.')[0][-1]]=(pickle.load(file(fname,'r')))
    memory={}    
    for key in DicList.keys():
        memory[key]=[]
        print key, 'iterated!'
        for sublist in DicList[key]:
            average1=0
            average2=0
            count=0 
            for Dict in sublist:
                print Dict
                for subkey, subval in Dict.iteritems():
                    if subval[0]!='UNIPROT':
                        if int(subval[0])<0:
                            print Dict
                        count+=1
                        average1+=int(subval[0])
                        average2+=int(subval[1])
            average1=float(average1)/float(count)
            average2=float(average2)/float(count)
            memory[key].append((average1,average2))
    srtd=sorted(memory.iteritems(), key=operator.itemgetter(0))
    for key, val in srtd:
        print key, val[0][0], val[1][0], val[2][0], '|', val[0][1], val[1][1], val[2][1]
            

def checkMatrix():
    pickleDump3=file('pickleDump3.dump','r')
    ValueMatrix=pickle.load(pickleDump3)
    NzeroList=ValueMatrix.nonzero()
    faultyList=[]
    for i in range(0,len(NzeroList[0])):
        index=(int(NzeroList[0][i]),int(NzeroList[1][i]))
        value=ValueMatrix[index[0],index[1]]
        if value < 0.0 or value > 1.0:
            faultyList.append(index)
    print len(faultyList)
    print len(NzeroList[0])
#     rsample=random.sample(range(0,len(NzeroList[0])), 10)
#     for i in rsample:
#         index=(int(NzeroList[0][i]),int(NzeroList[1][i]))
#         value=ValueMatrix[index[0],index[1]]
#         print index, value

def processEigenVectors():
    init=time()
    pickleDump=file('pickleDump.dump','r')
    eigenvals, eigenvects=pickle.load(pickleDump)
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    eigenVectsList=[]
    SuperIndex=[]
    print 'depicked',time()-init
    for i in range(0,eigenvects.shape[1]):
        eigenVectsList.append(eigenvects[:,i])
    i=0
    for eigenvect in eigenVectsList:
        i+=1
        normalized=np.multiply(eigenvect,eigenvect)
        indexes={}
        for k in range(0,10):
            index=np.argmax(normalized, 0)
            indexes[MatrixNumber2NodeID[index]]=(normalized[index], get_Descriptor(MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, index))
            normalized[index]=0
        SuperIndex.append(indexes)
    for subIndex in SuperIndex:
        for key in subIndex.keys():
            print '\t', key, subIndex[key]
        print '<=========================>'
    return SuperIndex


def columnSort():
    '''
    np.sum is broken for sparse matrixes. does the same thing with option axis=0
    '''
    pickleDump3=file('pickleDump3.dump','r')
    ValueMatrix=pickle.load(pickleDump3)
    SupportDict={}
    IndexDict={}
    nz=ValueMatrix.nonzero()
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    print 'ended imports'
    for i in range(0, len(nz[0])):
        if nz[0][i] not in SupportDict.keys():
            SupportDict[nz[0][i]]=[]
        SupportDict[nz[0][i]].append(nz[1][i])
    print 'ended loading indexes'
    for key in SupportDict.keys():
        if len(SupportDict[key])>1:
            IndexDict[key]=0
            for val in SupportDict[key]:
                IndexDict[key]+=float(ValueMatrix[key,val])
                
    srtd=sorted( IndexDict.iteritems(), key=operator.itemgetter(1),reverse=True )
    
    outf=file('columnSum.csv','w')
    
    for elt in srtd:
        Stri=str(str(MatrixNumber2NodeID[elt[0]])+'\t'+str(elt[1])+'\n')
        print  MatrixNumber2NodeID[elt[0]], elt[1]
        outf.write(Stri)
    outf.close()

# getMatrix(DfactorDict, 100)
#  
# checkMatrix()
# 
# get_eigenvect_Stats()

# processEigenVectors()

mass_Calibrate(6,10,True)
 
treat_Calibration()

# columnSort()

# TODO: create GO and Pathway Structure access
# Calibrate the values so that after ~ 3 transitions the correlation vanishes on average (Follow Pamela Silver Approach) => this is actually the cumulated perturbation of
# two targets that shoudl vanish totally 

# Shut down HiNT analysis => Slightly improves the result

# Synchronious eigenvectors approach: protect agains entering into a forbidden list the target node
# start iterating matrix multiplications starting from the node1 to go to the node2
# enter each node visited in the forbidden set, except for node2
# terminate iterating when there are no more new reaches for node2 after all the interations

# Percentage of information reaching a given node compared to all the information reaching the node: eigenvalue approach too.
# Error we do: compute three times

# Ok, what is going on is that we have collections of ~ 300 elements completely screwing our system

# The problem that a information broadcasting between the elements of the same group is not a good thing, but a direct broadcasting into a reaction is actually
# what we need in our matrix.