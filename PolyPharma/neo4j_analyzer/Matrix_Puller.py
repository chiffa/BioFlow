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
from scikits.sparse.cholmod import cholesky
from scipy.sparse import csc_matrix
from scipy.sparse import diags
from math import  sqrt
from Utils.Prot_Aboundances import ID2Aboundances
from scipy.sparse.csgraph import connected_components
from scipy.sparse.csgraph import shortest_path
import pylab

edge_type_filter0_1=["is_part_of_collection"]
edge_type_filter0_2=["is_same"]
edge_type_filter1_1=["is_Catalysant","is_reaction_particpant"]
edge_type_filter1_2=["is_part_of_complex", "is_Regulant"]
edge_type_filter1_3=["is_interacting"]
edge_type_filter2=["is_possibly_same"]

DfactorDict={"Group":0.5,
             "Same":1,
             "Reaction":0.33,
             "Contact_interaction":0.33,
             "possibly_same":0.1,
             }

ConductanceDict={"Group":20,
             "Same":100,
             "Reaction":1,
             "Contact_interaction":1,
             "possibly_same":0.1,
             }

# edge_type_filter3=["is_part_of_pathway","is_next_in_pathway"]

# Matrix filled with 0 (full dissipation, no connection), 1 (no dissipation, aliases) or r in [0,1] for partial dissipations



ReactionsList=[DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]

def get_Reaction_blocks(Connexity_Aware):
    ReagentClusters=[]
    Seeds=set()
    count=0
    for ReactionType in ReactionsList:
        for Reaction in ReactionType.get_all():
            Connex=True
            if Connexity_Aware:
                Connex=False
                if 'Main_Connex' in Reaction.custom:
                    Connex=True 
            if Connex and Reaction!=None:
                LocalList=[]
                for edge_type in edge_type_filter1_1:
                    if Reaction.bothV(edge_type)!=None:
                        for elt in Reaction.bothV(edge_type):
                            ID=str(elt).split('/')[-1][:-1]
                            if ID not in IDFilter:
                                LocalList.append(ID)
                                count+=1
                if len(LocalList)>1:
                    ReagentClusters.append(copy.copy(LocalList))
                    Seeds.update(LocalList)
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
        if len(LocalList)>0:
            Clusters[element]=copy.copy(LocalList)
    return Clusters, SuperSeed, count

def compute_Uniprot_Attachments():
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    Uniprot_Attach={} #Attaches the Uniprots to the proteins from the reactome, allowing to combine the information circulation values
    for SP_Node_ID in Uniprots:
        Vertex=DatabaseGraph.UNIPORT.get(SP_Node_ID)
        Generator=Vertex.bothV("is_same")
        if Generator!=None:
            Uniprot_Attach[SP_Node_ID]=[]
            for item in Generator:
                ID=str(item).split('/')[-1][:-1]
                Uniprot_Attach[SP_Node_ID].append(ID)
        print SP_Node_ID, 'processed', len(Uniprot_Attach[SP_Node_ID])
    pickle.dump(Uniprot_Attach,file('UP_Attach.dump','w'))
    return Uniprot_Attach

def load_Uniprot_Attachments():
    Uniprot_Attach=pickle.load(file('UP_Attach.dump','r'))
    return Uniprot_Attach    

def build_correspondances(IDSet,Rapid):
    NodeID2MatrixNumber={}
    MatrixNumber2NodeID={}
    ID2displayName={}
    ID2Type={}
    Uniprots=[]
    ID2Localization={}
    counter=0
    LocationBufferDict={}
    if not Rapid:
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
        pickleDump2=file('pickleDump2.dump','w')
        pickle.dump((NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots),pickleDump2)
        pickleDump2.close()
    else:
        NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
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


def getMatrix(decreaseFactorDict, numberEigvals, FastLoad, ConnexityAwareness=True):
    init=time()
    # Connect the groups of ingredients that share the same reactions1
    # Retrieve seeds for the matrix computation
    ReactLinks, InitSet,GroupLinks, GroupSet,SecLinks, SecSet,UP_Links, UPSet,HiNT_Links, FullSet,Super_Links,ExpSet=({},[],{},[],{},[],{},[],{},[],{},[])
    if not FastLoad:
        ReactLinks, InitSet, c = get_Reaction_blocks(ConnexityAwareness)
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
        
        for i in range(0,5):
            SecLinks2, SecSet2, c = get_expansion(SecSet,edge_type_filter1_2)
            print len(SecLinks2), len(SecSet2), c
            print time()-init, time()-t
            SecSet=SecSet2
            SecLinks=SecLinks2
            t=time()
        
        UP_Links, UPSet, c = get_expansion(SecSet,edge_type_filter0_2)
        print len(UP_Links), len(UPSet), c
        print time()-init, time()-t
#         t=time()
#         HiNT_Links, FullSet, c = get_expansion(UPSet,edge_type_filter1_3)
#         print len(HiNT_Links), len(FullSet), c
#         print time()-init, time()-t
#         t=time()
#         Super_Links, ExpSet, c = get_expansion(FullSet,edge_type_filter2)
#         print len(Super_Links), len(ExpSet), c
#         print time()-init, time()-t
        t=time()
        DF=file('dump5.dump','w')
#         pickle.dump((ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet),DF)
        pickle.dump((ReactLinks, InitSet,GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet),DF)
    else:
        DF=file('dump5.dump','r')
#         ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links,ExpSet=pickle.load(DF)
        ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet = pickle.load(DF)
    
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
    
    print "building correspondances", time()-init,
    t=time()
    # In_Case of the full impact: replace UPSet by ExpSet
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = build_correspondances(UPSet,FastLoad)
    
    
    print "building the ValMatrix", time()-init, time()-t
    t=time()
    
    # Group node definintion have to be corrected so they are not all related together but instead are linked towards the central "group" node!!!!
    
    # TODO: update all this element so that it incorporates the new linkage 
    
    # In_Case of the full impact: replace UPSet by ExpSet
    loadLen=len(UPSet)
    ValueMatrix=lil_matrix((loadLen,loadLen))
    ConductanceMatrix=lil_matrix((loadLen,loadLen))
    for group in ReactLinks:
        for elt in itertools.permutations(group,2):
            element=(NodeID2MatrixNumber[elt[0]],NodeID2MatrixNumber[elt[1]])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Reaction"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Reaction"]
            ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["Reaction"]/2.0
            ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["Reaction"]/2.0
            # the 2.0 is there because each symmetrical position will be attained twice 
            
    for key in GroupLinks.keys():
        for val in GroupLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Group"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Group"]
            
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Group"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Group"]
            
            ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["Group"]
            ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["Group"]
            
    for key in SecLinks.keys():
        for val in SecLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Contact_interaction"]
            
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Contact_interaction"]
            
            ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["Contact_interaction"]
            ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["Contact_interaction"]
            
    for key in UP_Links.keys():
        for val in UP_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Same"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Same"]
            
            element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
            ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Same"],1)
            ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Same"]
            
            ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["Same"]
            ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["Same"]
            
#     for key in HiNT_Links.keys():
#         for val in HiNT_Links[key]:
#             element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
#             ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
#             ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Contact_interaction"]
#             
#             element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
#             ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["Contact_interaction"],1)
#             ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["Contact_interaction"]
#             
#             ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["Contact_interaction"]
#             ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["Contact_interaction"]
#             
#     for key in Super_Links.keys():
#         for val in Super_Links[key]:
#             element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
#             ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["possibly_same"],1)
#             ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["possibly_same"]
#             
#             element=(NodeID2MatrixNumber[val],NodeID2MatrixNumber[key])
#             ValueMatrix[element[0],element[1]]=min(ValueMatrix[element[0],element[1]]+decreaseFactorDict["possibly_same"],1)
#             ConductanceMatrix[element[0],element[1]]=ConductanceMatrix[element[0],element[1]]-decreaseFactorDict["possibly_same"]
#             
#             ConductanceMatrix[element[1],element[1]]=ConductanceMatrix[element[1],element[1]]+decreaseFactorDict["possibly_same"]
#             ConductanceMatrix[element[0],element[0]]=ConductanceMatrix[element[0],element[0]]+decreaseFactorDict["possibly_same"]

    print "entering eigenvect computation", time()-init, time()-t
    t=time()
    eigenvals, eigenvects = eigsh(ValueMatrix,numberEigvals)
    eigenvals2, eigenvects2 = eigsh(ConductanceMatrix,numberEigvals)
    print eigenvals
    print '<======================>'
    print eigenvals2
    print '<======================>'
    print np.all(eigsh(ConductanceMatrix)[0] > 0)
    output = file('eigenvals.csv','w')
    output.write(eigenvals)
    output.close()
    output2 = file('eigenvals.csv','w')
    output2.write(eigenvals2)
    output2.close()
    pickleDump=file('pickleDump.dump','w')
    pickle.dump((eigenvals, eigenvects), pickleDump)
    pickleDump.close()
    pickleDump_0_5=file('pickleDump_0_5.dump','w')
    pickle.dump((eigenvals2,eigenvects2),pickleDump_0_5)
    pickleDump_0_5.close()
    pickleDump3=file('pickleDump3.dump','w')
    pickle.dump(ValueMatrix,pickleDump3)
    pickleDump3.close()
    pickleDump4=file('pickleDump4.dump','w')
    pickle.dump(ConductanceMatrix,pickleDump4)
    pickleDump4.close()
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

def get_voltages(numpy_array, MatrixNumber2NodeID, InformativityDict):
    for i in range(0, len(numpy_array)):
        InformativityDict[MatrixNumber2NodeID[i]]+=numpy_array[i,0]
    return
        
def create_InfoDict(MatrixNumber2NodeID):
    new_dict={}
    for val in MatrixNumber2NodeID.values():
        new_dict[val]=0.0
    return new_dict

def compute_sample_circulation_intensity(Sample, epsilon=1e-10, array_v=''):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    @param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    @type: list of ints
    
    @param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    @type epsilon: float
    
    '''
    init=time()
    conductance_Matrix=pickle.load(file('pickleDump4.dump','r'))
    # NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    InformativityArray=np.zeros((conductance_Matrix.shape[0],1))                                # Database ID to Informativity
    Solver=cholesky(csc_matrix(conductance_Matrix),epsilon)
    # The informativities are calculated only by using the uniprot proteins as the source and extraction points.
    # This reduces the number of interations from  25k to 5 and the number of LU decompositions in a similar manner
    print 1, time()-init
    init=time()
    print Sample
    print len(Sample)
    for i in range(0,len(Sample)):
        print '\n', 2, i, time()-init
        init=time()
        for j in range(i,len(Sample)):
            if j%25==24:
                print '*',
            J=np.zeros((conductance_Matrix.shape[0],1))
            J[Sample[i],0]=1.0
            J[Sample[j],0]=-1.0
            V=Solver(J)
            Current=get_Current_all(conductance_Matrix,V,J)
            InformativityArray+=Current
            name='InfoArray_new'+str(array_v)+'.dump'
    pickle.dump(InformativityArray,file(name,'w'))
    return InformativityArray

def Compute_circulation_intensity(epsilon=1e-10):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    @param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    @type epsilon: float
    
    '''
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    re_UP=[]
    for SP in Uniprots:
        re_UP.append(NodeID2MatrixNumber[SP])
    return compute_sample_circulation_intensity(re_UP,epsilon,'_Full')


def Compute_random_sample(sample_size, iterations, epsilon=1e-10, name_version=''):
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    re_UP=[]
    for SP in Uniprots:
        re_UP.append(NodeID2MatrixNumber[SP])
    for i in range(0,iterations):
        random.shuffle(re_UP)
        re_UP=re_UP[:sample_size]
        compute_sample_circulation_intensity(re_UP,epsilon,name_version+str(i))
    
    
def Analyze_relations(Current, number,MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization,AdditionalInfos=''):
    '''
    Just a way to print out the nodes making pass the most of the current
    
    @param Current: the current
    @type Current: numpy array
      
    @param number: the number of the most important relations one wants to observe
    @type number: integer
    
    '''
    writer=file('Current_Analysis.csv','w')
    infos=np.absolute(Current)
    initLine='ID'+'\t'+'Mean Informativity'+'\t'+'Type'+'\t'+'Standard deviation'+'\t'+'Pessimistic Info estimation'+'\t'+'optimistic Info estimation'+'\t'+'Protein Aboundace'+'\t'+'Essential for Overington'+'\t'+'displayName'+'\t'+'localization'+'\n'
    writer.write(initLine)
    if AdditionalInfos=='':
        for i in range(0,min(number,len(infos))):
            locmax=np.argmax(infos)
            Buffer=''
            Buffer+=str(MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
#             print MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer+=str(ID2Type[MatrixNumber2NodeID[locmax]])+'\t'
#             print ID2Type[MatrixNumber2NodeID[locmax]], '\t',
            Buffer+=str(ID2displayName[MatrixNumber2NodeID[locmax]])+'\n'
#             print ID2displayName[MatrixNumber2NodeID[locmax]]
            infos[locmax]=0.0
            writer.write(Buffer)
    else:
        for i in range(0,min(number,len(infos))):
            locmax=np.argmax(infos)
            Buffer=''
            Buffer+=str(MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
#             print MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer+=str(ID2Type[MatrixNumber2NodeID[locmax]])+'\t'
#             print ID2Type[MatrixNumber2NodeID[locmax]], '\t',
            for elt in AdditionalInfos[locmax,:].tolist():
                Buffer+=str(elt)+'\t'
#                 print elt, '\t',
            Buffer+=str(ID2displayName[MatrixNumber2NodeID[locmax]])
            if MatrixNumber2NodeID[locmax] in ID2Localization.keys():
                Buffer+='\t'+str(ID2Localization[MatrixNumber2NodeID[locmax]])
            Buffer+='\n'
#             print ID2displayName[MatrixNumber2NodeID[locmax]]
            infos[locmax]=0.0
            writer.write(Buffer)
    writer.close()

def Info_circulation_for_Single_Node(conductance_Matrix,source_MatrixID,sinks_MatrixIDs,epsilon=1e-10):
    '''
    Computes the information circulation for a single node, according to a slightly improved method compared to the one described by Missiuro
    @param conductance_Matrix: Conductance matrix
    @type conductance_Matrix: scipy.sparse lil_matrix
    
    @param source_MatrixID: number of the row/line of the source node within the conductance matrix corresponding to the source
    @type source_MatrixID: int
    
    @param sinks_MatrixIDs: list of numbers of the rows/lines within the conductance matrix corresponding to sinks  
    @type sinks_MatrixIDs: list of ints
    
    @param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    @type epsilon: float
    '''
    Solver=cholesky(csc_matrix(conductance_Matrix),epsilon)
    Cummulative_Informativity=np.zeros((conductance_Matrix.shape[0],1))
    for sink in sinks_MatrixIDs:
        J=np.zeros((conductance_Matrix.shape[0],1))
        J[source_MatrixID,0]=1.0
        J[sink,0]=-1.0
        Voltage=Solver(J)
        current=get_Current_all(conductance_Matrix,Voltage,J)
        Cummulative_Informativity+=current
    return Cummulative_Informativity

def get_Current_SingleNode(Matrix_NodeID,Conductance_Matrix,Voltages):
    '''
    Recovers the current within one single node
    
    @param Matrix_NodeID: number of the row/line of the source node within the conductance matrix corresponding to the source
    @type Matrix_NodeID: int
    
    @param Conductance_Matrix: Conductance matrix
    @type Conducntace_Matrix: scipy.sparse lil_matrix
    
    @param Voltages: Informativity gradient in each node obtained by the solution of the matrix equation ConductanceMatrix*Voltage = J
    @type Voltages: numpy array

    '''
    submatrix=Conductance_Matrix[:,Matrix_NodeID].todense()
    submatrix[Matrix_NodeID,0]=0.0
    basevoltage=np.ones(submatrix.shape)
    basevoltage=Voltages[Matrix_NodeID,0]*basevoltage
    Current=float(np.sum(np.absolute(np.multiply(basevoltage,submatrix)-np.multiply(Voltages,submatrix)))/2.0)
    return Current

def get_Current_NodeSet(NodeSet,Conductance_Matrix,Voltages):
    '''
    Recovers the current within a set of Nodes
    
    @param NodeSet: number of the row/line of the source node within the conductance matrix corresponding to the source
    @type NodeSet: list of ints
    
    @param Conductance_Matrix: Conductance matrix
    @type Conducntace_Matrix: scipy.sparse lil_matrix
    
    @param Voltages: Informativity gradient in each node obtained by the solution of the matrix equation ConductanceMatrix*Voltage = J
    @type Voltages: numpy array

    '''    
    
    Dico={}
    for elt in NodeSet:
        Dico[elt]=get_Current_SingleNode(elt,Conductance_Matrix,Voltages)
    return Dico

def get_Current_all(Conductance_Matrix, Voltages, J):
    '''
    Recovers the current for all the nodes
    
    @param Conductance_Matrix: Conductance matrix
    @type Conducntace_Matrix: scipy.sparse lil_matrix
    
    @param Voltages: Informativity gradient in each node obtained by the solution of the matrix equation ConductanceMatrix*Voltage = J
    @type Voltages: numpy array
    
    '''
    diag_Voltages=lil_matrix(diags(Voltages.T.tolist()[0],0))
    Corr_Conductance_Matrix=Conductance_Matrix-lil_matrix(diags(Conductance_Matrix.diagonal(),0))
    Currents=(np.absolute(diag_Voltages*Corr_Conductance_Matrix-Corr_Conductance_Matrix*diag_Voltages).sum(axis=0).T+np.absolute(J))/2.0
    return Currents

    
def stats_over_random_info_circ_samples(UniProtAttachement=True):
    UPNode_IDs_2Proteins_IDs_List=load_Uniprot_Attachments()
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    DicList=[]
    FnameList=[]
    Filenames = listdir('.')
    for filename in Filenames:
        if 'InfoArray_new' in filename:
            FnameList.append(filename)
    for fname in FnameList:
        DicList.append(pickle.load(file(fname,'r')))
    Stats_Mat=np.concatenate(tuple(DicList),axis=1)
    MeanInfos=np.mean(Stats_Mat,axis=1).reshape((Stats_Mat.shape[0],1))
    STDInfos=np.std(Stats_Mat,axis=1).reshape((Stats_Mat.shape[0],1))
    if UniProtAttachement:
        for UP_ID in Uniprots:
            StdBuffer=float(STDInfos[NodeID2MatrixNumber[UP_ID]])**2
            for Prot_ID in UPNode_IDs_2Proteins_IDs_List[UP_ID]:
                MeanInfos[NodeID2MatrixNumber[UP_ID],0]+=MeanInfos[NodeID2MatrixNumber[Prot_ID],0]
                StdBuffer+=float(STDInfos[NodeID2MatrixNumber[Prot_ID],0])**2
                ID2Localization[UP_ID]=ID2Localization[Prot_ID]
            STDInfos[NodeID2MatrixNumber[UP_ID],0]=sqrt(StdBuffer)
    pessimist=MeanInfos-1.97*STDInfos
    optimist=MeanInfos+1.97*STDInfos
    aboundances=np.zeros((Stats_Mat.shape[0],1))
    errcount1=0
    for key,val in ID2Aboundances.iteritems():
        if key in NodeID2MatrixNumber.keys():
            aboundances[NodeID2MatrixNumber[key],0]=val
        else:
            errcount1+=1
    Overingtonicitiy=np.zeros((Stats_Mat.shape[0],1))
    Overington_IDList=pickle.load(file('IDList.dump','r'))
    errcount2=0
    for ID in Overington_IDList:
        if ID in NodeID2MatrixNumber.keys():
            Overingtonicitiy[NodeID2MatrixNumber[ID],0]=1.0
        else:
            errcount2+=1
    print 'errcount from stats', errcount1, errcount2
    Additional=np.concatenate((STDInfos,pessimist,optimist,aboundances,Overingtonicitiy),axis=1)
    Analyze_relations(MeanInfos, 50000, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Additional)

def check_Silverality(sample_size,iterations):
    '''
    Checks how rapidly the information passing through the neighbouring proteins decreases
    with the distance between them
    '''
    from scipy.sparse.csgraph import dijkstra
    conductance_Matrix=pickle.load(file('pickleDump4.dump','r'))
    value_Matrix=pickle.load(file('pickleDump3.dump','r'))
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    cumulative=[]
    rei_UP=[]
    for SP in Uniprots:
        rei_UP.append(NodeID2MatrixNumber[SP])
    for i in range(0,iterations):
        random.shuffle(rei_UP)
        re_UP=rei_UP[:sample_size+1]
        source_MatrixID=re_UP[0]
        sinks_MatrixIDs=re_UP[0:]
        Currents=Info_circulation_for_Single_Node(conductance_Matrix,source_MatrixID,sinks_MatrixIDs,epsilon=1e-10)
        random.shuffle(rei_UP)
        targets=rei_UP[:sample_size/10]
        distances = dijkstra(value_Matrix, indices=source_MatrixID, unweighted=True)
        for index in targets:
            cumulative.append((distances[index],Currents[index,0]))
        for i in range(0,sample_size/5):
            index=np.argmin(distances)
            print distances
            cumulative.append((distances[index],Currents[index,0]))
            distances[index]=100
    pickle.dump(cumulative,file('Silverality.dump','w'))
    return cumulative

def analyze_Silverality():
    Silverality_List=pickle.load(file('Silverality.dump','r'))
    SuperDict={}
    for distance, value in Silverality_List:
        if distance not in SuperDict.keys():
            SuperDict[distance]=[]
        SuperDict[distance].append(value)
    
    for key, val in SuperDict.iteritems():
        print key, val
    
    ArrayDict={}
    StatsDict={}
    x=np.zeros((10,1))
    y1=np.zeros((10,1))
    y2=np.zeros((10,1))
    y3=np.zeros((10,1))
    i=0
    Srtd=sorted(SuperDict.iteritems(),key=operator.itemgetter(0))
    for key, val in Srtd:
        if i>9:
            break
        ArrayDict[key]=np.array(val)
        StatsDict[key]=(np.mean(ArrayDict[key]),np.std(ArrayDict[key]))
        x[i,0]=key
        y1[i,0]=StatsDict[key][0]
        y2[i,0]=StatsDict[key][0]+StatsDict[key][1]
        y3[i,0]=StatsDict[key][0]-StatsDict[key][1]
        i+=1

    
    pylab.plot(x, y1, '-k', label='mean')
    pylab.plot(x, y2, '-r', label='mean+std')
    pylab.plot(x, y3, '-b', label='mean-std')
    pylab.legend(loc='upper right')
    pylab.show()

#     array1=np.zeros((len(Silverality_List),1))
#     array2=np.zeros((len(Silverality_List),1))
#     
#     for i in range(0, len(Silverality_List)):
#         array1[i,0]=Silverality_List[i][0]
#         array2[i,0]=Silverality_List[i][1]
#     pylab.plot(array1, array2)
#     pylab.show()
    # TODO: perform statistical analysis

def Perform_Loading_Routines():
    '''
    @param eigenvals:
    @type eigenvals:
    
    @param param:   
    '''
    
    getMatrix(DfactorDict, 100, False, False,)
    compute_Uniprot_Attachments()
    
    Write_Connexity_Infos()

def Perform_Testing_Routines():
    raise NotImplementedError

def Write_Connexity_Infos():    
    ValueMatrix = pickle.load(file('pickleDump3.dump','r'))
    CCOmps = connected_components(ValueMatrix, directed=False)
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    
    counters = np.zeros((CCOmps[0],1))
    for elt in range(0,len(CCOmps[1])):
        counters[CCOmps[1][elt],0] += 1
    major_Index=np.argmax(counters)
    
    ln=len(CCOmps[1])
    for i in range(0,len(CCOmps[1])):
        print 'running,',i,ln
        if CCOmps[1][i]==major_Index:
            Node=DatabaseGraph.vertices.get(MatrixNumber2NodeID[i])
            Node.custom='Main_Connex'
            Node.save()

def Erase_Additional_Infos():
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    
    for NodeID in NodeID2MatrixNumber.keys():
        Node=DatabaseGraph.vertices.get(NodeID)
        Node.custom=''
        Node.save()


# TODO: reverse GO_Access: provided the Uniprots find the proteins carrying over the most information

stats_over_random_info_circ_samples(True)
# 
# check_Silverality(100,100)
# 
# analyze_Silverality()
# 
# Write_Connexity_Infos()

# getMatrix(DfactorDict, 100, False)

# compute_Uniprot_Attachments()

# Compute_circulation_intensity()

# Compute_random_sample(100,3,1e-10, '3_')

# stats_over_random_info_circ_samples(True)

# this is going to last for a while. Now we need to get it split among several nodes

# Compute_circulation_intensity()
#  
# checkMatrix()
# 
# get_eigenvect_Stats()

# processEigenVectors()

# mass_Calibrate(6,10,True)
#  
# treat_Calibration()

# columnSort()



# TODO: There might be an error in the module responsible for linkage between the uniprots and the accession numbers: for instance the 20253 has an annotation with an Acnum, but
# has no Uniprot attached to it within the database

# TODO: create GO and Pathway Structure access
# Calibrate the values so that after ~ 3 transitions the correlation vanishes on average (Follow Pamela Silver Approach) => this is actually the cumulated perturbation of
# two targets that shoudl vanish totally 

# TODO: implementation while using Ehit interactions only


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


# In order to be precise, we should not only take in account the power of bindinb between a molecule and protein and criticality of the protein, but also the abundance of the
# protein in the reactome

# => Done with the aboundance retrieval

# DONE: use sparse matrixes routines to calculate the number of connex elements in the graph
#   Problem: there are 58 disconnected sets.
#   Solution: retrieve the Node Ids of the main connex Set and write them into the neo4j graph, then retrieve only them

# DONE: markup of the major connex graph within neo4j database
#    Waiting for the execution


# DONE: calculate the distance graph
    # seems to work pretty well with Djikistra.
    # Can we perform a retrieval of specific nodes within distance X of the main component? 

# DONE: buid jump tables to compute the number of reactional transitions
#    Implemented by using djikstra algo from scipy.sparse.csgraph
#

# DONE: retrieve Pamela silver's degradation of the data with the time
#    Waiting for the execution
#    

# DONE: pull in the annotations regarding the proteins aboundances
#    
#    

# DONE: pull in the 300 essential targets from the EBI dude (John Overington)
#     Results aren't so conclusive. It seems that the protein concentration defenitely plays some role in the determining if a protein is a 
#     Target of an existing drug or not, butthe informativity seems not. Probably this is due to the fact that the targeted proteins are often 
#     cellular receptors.

# DONE: perform a localization factor pull-out for the Uniprots based on their proteins of attachement
#        Waiting for the execution

# DONE: broadcast to uniprots for the localization of the pointed proteins