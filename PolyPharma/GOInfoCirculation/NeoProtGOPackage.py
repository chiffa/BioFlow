'''
Created on Mar 27, 2013

@author: akucahravy
'''

import os
os.environ['NEO4J_PYTHON_JVMARGS']='-Xms128M -Xmx512M'
os.environ['JAVA_HOME']='/usr/lib/jvm/java'
from neo4j import GraphDatabase

from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import func, select
from time import time
from TargetPreProcessing import HumanTargets_proc
import Useful_Dicts as UD
import datetime
import math
import operator
import itertools
import numpy as np
import random as rand

np.set_printoptions(precision=2)

lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb', echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

 
neo=GraphDatabase('/home/akucahravy/DB/neo4j')


#with neo.transaction:
#    proteins=neo.node()
#    gos=neo.node()
#    neo.reference_node.PROTEINS(proteins)
#    neo.reference_node.GOS(gos)
#    prot_idx=neo.node.indexes.create('proteins')
#    gos_idx=neo.node.indexes.create('gos')
    
with neo.transaction:
    prot_idx=neo.node.indexes.get('proteins')
    gos_idx=neo.node.indexes.get('gos')
    for rel in neo.reference_node.GOS:
        gos = rel.end
    for rel in neo.reference_node.PROTEINS:
        proteins = rel.end

ProtNamesBuffer=set()
GONamesBuffer=set()
EdgesBuffer=set()

def create_Protein(name,Importance):
    if prot_idx['name'][name].single!=None:
        return prot_idx['name'][name].single
    else:
        with neo.transaction:
            ProtObj = neo.node(name=name,imp=Importance)
            ProtObj.INSTANCE_OF(proteins)
            prot_idx['name'][name]=ProtObj 
        return ProtObj

def create_GO(GOId, Informativity):
    if gos_idx['name'][GOId].single!=None:
        return gos_idx['name'][GOId].single
    else:
        with neo.transaction:
            GOObj = neo.node(name=GOId,inf=Informativity)
            GOObj.INSTANCE_OF(gos)
            gos_idx['name'][GOId]=GOObj 
        return GOObj

def link(ProtName,GOName):
    GOObj=get_GOObj(GOName)
    ProtObj=get_ProtObj(ProtName)
    Duplicate=False
    for relationship in ProtObj.DESCRIBES:
        if relationship.start == GOObj:
            Duplicate=True
    if not Duplicate:
        with neo.transaction:
            GOObj.DESCRIBES(ProtObj)
    return GOObj

def batchCreate_links(EdgeList):
    # EdgeList is of the form {(ProtName,GOName):(imp,inf)}
    Prots_TBC={}
    GOs_TBC={}
    Edges_TBC=[] # edges to be created
    Edges_TC=[] # edges to check
    for key in EdgeList.keys():
        if key[0] not in ProtNamesBuffer:
            if prot_idx['name'][key[0]].single==None:
                Prots_TBC[key[0]]=EdgeList[key][0]
                Edges_TBC.append(key)
            else:
                Edges_TC.append(key)
                ProtNamesBuffer.add(key[0])
        else:
            Edges_TC.append(key)
        if key[1] not in GONamesBuffer:
            if gos_idx['name'][key[1]].single==None:
                GOs_TBC[key[1]]=EdgeList[key][1]
                Edges_TBC.append(key)
            else:
                Edges_TC.append(key) 
                GONamesBuffer.add(key[1])
        else:
            Edges_TC.append(key)
            
    Edges_TC=set(Edges_TC)
    Edges_TBC=set(Edges_TBC)
    Edges_TC=Edges_TC-Edges_TBC
    Edges_TC=Edges_TC-EdgesBuffer
    Edges_TC=list(Edges_TC)
    Edges_TBC=list(Edges_TBC)
    
    #Control the Edges_TC:
    
    # Group according to the ProtObj,i.e. first value
    ProtNames=[]
    for entry in Edges_TC:
        ProtNames.append(entry[0])
    
    ProtNames=list(set(ProtNames))
    
    for ProtName in ProtNames:
        sublist=[]
        for key in Edges_TC:
            if key[0]==ProtName:
                sublist.append(key)
        ProtObj=get_ProtObj(ProtName)
        
        GoIds=[]
        for key in sublist:
            GoIds.append(key[1])
        
        for relationship in ProtObj.DESCRIBES:
            nom=relationship.start['name']
            if nom in GoIds:
                EdgesBuffer.add((ProtName,nom))
    
    Edges_TC=list(set(Edges_TC)-EdgesBuffer)
    Edges_TBC=list(set(Edges_TC)|set(Edges_TBC))
    
    #Here goes the transaction!!!!
    with neo.transaction:
        #Create all the Nodes:
        for ProtName, importance in Prots_TBC.items():
            ProtObj = neo.node(name=ProtName,imp=importance)
            ProtObj.INSTANCE_OF(proteins)
            prot_idx['name'][ProtName]=ProtObj
        for GOId, informativity in GOs_TBC.items():
            GOObj = neo.node(name=GOId,inf=informativity)
            GOObj.INSTANCE_OF(gos)
            gos_idx['name'][GOId]=GOObj
        #Create all the edges:
        for ProtName,GoName in Edges_TBC:
            GOObj=get_GOObj(GoName)
            ProtObj=get_ProtObj(ProtName)
            GOObj.DESCRIBES(ProtObj)       
    

def get_ProtObj(ProtName):
    return prot_idx['name'][ProtName].single

def get_GOObj(GOId):
    return gos_idx['name'][GOId].single

def copy_N_entries(number):
    Query=session.query(TableBuilder.UNIPROT_FullGO2).\
                            limit(number).\
                            all()
    Prots=[]
    GOs=[]
    links=[]
    for row in Query:
        Prots.append(str(row.UNIPROTID))
        GOs.append(str(row.GOs))
        links.append((str(row.UNIPROTID),str(row.GOs))) 
    Prots=list(set(Prots))
    GOs=list(set(GOs))
    links=list(set(links))
    for ProtName in Prots:
        create_Protein(ProtName, 1.0)
    for GOId in GOs:
        create_GO(GOId, 1.0)
    for Llink in links:
        link(Llink[0],Llink[1])        
 

def copy_from_to_entries(offset,chunk):
    Query=session.query(TableBuilder.UNIPROT_FullGO2).\
                            offset(offset*chunk).\
                            limit(chunk).\
                            all()
    links={}
    for row in Query:
        links[(str(row.UNIPROTID),str(row.GOs))]=(1.0,1.0)
    batchCreate_links(links)

def load_all_chunks(chunk_size):
    Query=session.query(TableBuilder.UNIPROT_FullGO2).count()
    number_of_chunks=Query/chunk_size
    print 'total chunks:', number_of_chunks
    estimatedTime=100
    for i in range(0,number_of_chunks+1):
        start=time()
        copy_from_to_entries(i,chunk_size)
        estimatedTime=0.75*estimatedTime+0.25*(time()-start)
        print 'chunk', i, 'finished in', "{0:.2f}".format(time()-start), 'estimated chunks per second:',"{0:.2f}".format(chunk_size/estimatedTime) ,'Total Runtime',datetime.timedelta(seconds=estimatedTime*number_of_chunks)
        print '=================================================================>>>>>'

def set_batch_Prot_importance(ImportanceDict):
    # The argument is expected to be of the form {'ProtName','importance'}
    with neo.transaction:
        for ProtName, Importance in ImportanceDict.items():
            ProtObj=prot_idx['name'][ProtName].single
            ProtObj['imp']=Importance
  
def set_batch_GO_informativity(InfoDict):
    # The argument is expected to be of the form {'GOId','importance'}
    with neo.transaction:
        for GOId, informativity in InfoDict.items():
            GOObj=gos_idx['name'][GOId].single
            GOObj['inf']=informativity

def set_GO_inf_all(chunkSize):
    TotalChunks=session.query(TableBuilder.GO_Informativity).\
                        count()/chunkSize
    start=100
    print 'total Chunks:', TotalChunks                   
    for i in range(0,TotalChunks):
        print 'loading Chunk', i, 'estimated chunks per second:', "{0:.2f}".format(chunkSize/(time()-start))
        CurrentChunk=session.query(TableBuilder.GO_Informativity).\
                                offset(i*chunkSize).\
                                limit(chunkSize)
        LoadDict={}
        
        for row in CurrentChunk:
            LoadDict[row.go]=math.log(row.informativity,2)
        
        set_batch_GO_informativity(LoadDict)
        start=time()

def set_GO_names_for_all(chunkSize):
    Query=session.query(TableBuilder.UNIPROT_FullGO2.GOs).\
                        distinct().\
                        all()
    Valid_GOs=[]
    for elt in Query:
        Valid_GOs.append(elt[0])
    
    TotalChunks=session.query(TableBuilder.GO_Term).\
                        count()/chunkSize
    start=100
    print 'total Chunks:', TotalChunks                   
    for i in range(0,TotalChunks):
        print 'loading Chunk', i, 'estimated elts per second:', "{0:.2f}".format(chunkSize/(time()-start))
        CurrentChunk=session.query(TableBuilder.GO_Term).\
                                offset(i*chunkSize).\
                                limit(chunkSize)

        with neo.transaction:
            for row in CurrentChunk:
                if row.termid in Valid_GOs:
                    GOObj=gos_idx['name'][str(row.termid)].single
                    GOObj['fullname']=str(row.name)
                    GOObj['namespace']=str(row.namespace)
        start=time()


def export_all_to_gdf():
    print 'Oh-oh, export_all_to_gdf() does nothing yet'
    # send the gdf data to the 

def naive_pairwise_info_circ(ProtName1,ProtName2):
    #Batching it may increase performance later
    ProtObj1=prot_idx['name'][ProtName1].single
    ProtObj2=prot_idx['name'][ProtName2].single
    GOList1=[]
    GOList2=[]
    
    

    for relationship in ProtObj1.DESCRIBES:
        GOList1.append(relationship.start)
    
    for relationship in ProtObj2.DESCRIBES:
        GOList2.append(relationship.start)
    
    CommonGOs=list(set(GOList1) & set(GOList2))
    print 'CommonGOs are:'
    for GO in CommonGOs:
        print GO['name'], GO['inf']
    
    inf=0.0
    for GOObj in CommonGOs:
        inf+=GOObj['inf']
    
    print '<==================================>'
    return inf*(ProtObj1['imp']+ProtObj2['imp'])

def naive_percetaged_info_circ(ProtName1,ProtName2):
    #Batching it may increase performance later
    #Does 
    ProtObj1=prot_idx['name'][ProtName1].single
    ProtObj2=prot_idx['name'][ProtName2].single
    GOList1=[]
    GOList2=[]

    for relationship in ProtObj1.DESCRIBES:
        GOList1.append(relationship.start)
    
    for relationship in ProtObj2.DESCRIBES:
        GOList2.append(relationship.start)
    
    CommonGOs=list(set(GOList1) & set(GOList2))
    print 'CommonGOs are:'
    for GO in CommonGOs:
        print GO['name'], GO['inf']
    
    inf=0.0
    for GOObj in CommonGOs:
        inf+=GOObj['inf']
    
    print '<==================================>'
    return inf*(ProtObj1['imp']+ProtObj2['imp'])



def simple_rapid_naive_group_info(ProtNames2Importance):
    print 'starting!'
    ProtName2Obj={}
    ProtName2GOObj={}
    GOName2Importance={}
    Importance_Distribution={}
    for ProtName in ProtNames2Importance.keys():
        ProtName2Obj[ProtName]=prot_idx['name'][ProtName].single
        ProtName2GOObj[ProtName]=[]
        for relationship in ProtName2Obj[ProtName].DESCRIBES:
            ProtName2GOObj[ProtName].append(relationship.start)
    for pair in itertools.combinations(ProtName2Obj.keys(),2):
        # We are not doing self-alignement for each protein present
        CommonGOs=list(set(ProtName2GOObj[pair[0]]) & set(ProtName2GOObj[pair[1]]))
        #Check if the GO object is immutable and can be used as a key in a dictionary 
        for GO in CommonGOs:
            #Try to unbuffer this part and see if hte performance goes down drastically
            if not GO in GOName2Importance.keys():
                GOName2Importance[GO]=0.0
                Importance_Distribution[GO]={}
            Increase=(ProtNames2Importance[pair[0]]+ProtNames2Importance[pair[1]])*GO['inf']

            GOName2Importance[GO]+=Increase
            Importance_Distribution[GO][pair]=Increase
    
    sortedDict=sorted(GOName2Importance.iteritems(),key=operator.itemgetter(1),reverse=True)
    
    for GOObj, flow in sortedDict:
        if flow>0.01:   
            SortedSubdict=sorted(Importance_Distribution[GOObj].iteritems(),key=operator.itemgetter(1),reverse=True)
            Buffer=''
            for Pair, value in SortedSubdict:
                Buffer+=' |'+str(Pair[0])+'<->'+str(Pair[1])+': '+str("{0:.2f}".format(value/flow*100))+'%'
            print GOObj['name'], "{0:.2f}".format(GOObj['inf']), "{0:.2f}".format(flow), GOObj['namespace'], GOObj['fullname']
            print Buffer
            print '<==============================>'

# TODO: eliminate importance from the protein description in the neo4j database        
# TODO: Add a distinct method for inter-group info circulation

def rapid_naive_group_pairwise_info(ProtNames2Importance):
    # Used more to break into essential compounds
    # Should we cast the self-informativity to 0 or not? => Yes
    size=len(ProtNames2Importance)
    outmatrix=np.zeros((size,size))
    
    ProtName2Obj={}
    ProtName2GOObj={}
    Prot2Prot_Flow={}
    ProtName2Index={}
    ProtIndex=[]
    i=0
    for ProtName in ProtNames2Importance.keys():
        ProtName2Index[ProtName]=i
        ProtIndex.append(ProtName)
        i+=1
        ProtName2Obj[ProtName]=prot_idx['name'][ProtName].single
        ProtName2GOObj[ProtName]=[]
        for relationship in ProtName2Obj[ProtName].DESCRIBES:
            ProtName2GOObj[ProtName].append(relationship.start)
    for pair in itertools.combinations(ProtName2Obj.keys(),2): 
        # We are not doing self-alignement for each protein present
        CommonGOs=list(set(ProtName2GOObj[pair[0]]) & set(ProtName2GOObj[pair[1]]))
        #Check if the GO object is immutable and can be used as a key in a dictionary 
        if not pair in Prot2Prot_Flow.keys():
            Prot2Prot_Flow[pair]=0.0
        
        for GO in CommonGOs:
            #Try to unbuffer this part and see if hte performance goes down drastically
            Prot2Prot_Flow[pair]+=(ProtNames2Importance[pair[0]]+ProtNames2Importance[pair[1]])*GO['inf']
        outmatrix[ProtName2Index[pair[0]]][ProtName2Index[pair[1]]]=Prot2Prot_Flow[pair]
        outmatrix[ProtName2Index[pair[1]]][ProtName2Index[pair[0]]]=Prot2Prot_Flow[pair]
    
    return [ProtIndex,outmatrix]

def rapid_naive_group_pairwise_info_normalized(ProtNames2Importance):
    # Used more to break into essential compounds
    # This version performs the notmalization regarding how much similarity flow there are '
    # Between the proteins compared to the on that have the  ()
    
    ## TODO: check that it actually works as intended and that the
    #  Information flow is not all that skrewed all
    #  We want self-similarity to be equal to 100% and everyting else build from it.
    
    ## TODO: perform a more violent buffering of GOname2informativities to avoid redundant
    #  Calls to graph db
    
    # Should we cast the self-informativity to 0 or not? => Yes
    
    TotalInfoQuery=session.query(TableBuilder.UNIPROT_SpecInf).\
                            filter(TableBuilder.UNIPROT_SpecInf.TotalInfo<20.0).\
                            filter(TableBuilder.UNIPROT_SpecInf.UNIPROTID.in_(ProtNames2Importance.keys())).\
                            all()
    
    ProtNames_insuff_Info=[]
    
    for row in TotalInfoQuery:
        ProtNames_insuff_Info.append(str(row.UNIPROTID))
        del ProtNames2Importance[str(row.UNIPROTID)]
    
    print 'Following proteins do not follow Zipf distribution and were eliminated from analysis: ', len(ProtNames_insuff_Info), ProtNames_insuff_Info
    
    TotalInfoQuery=session.query(TableBuilder.UNIPROT_SpecInf).\
                            filter(TableBuilder.UNIPROT_SpecInf.UNIPROTID.in_(ProtNames2Importance.keys())).\
                            all()
    
    ProtNames_Yes_Info=[]
    
    for row in TotalInfoQuery:
        ProtNames_Yes_Info.append(str(row.UNIPROTID))
    
    No_set=list(set(ProtNames2Importance.keys())-set(ProtNames_Yes_Info))
    
    for elt in No_set:
        del ProtNames2Importance[elt]
    
    print 'Following proteins are not in SwissProt: ', len(No_set), No_set
    
    print 'Remaining proteins:', len(ProtNames2Importance.keys()), ProtNames2Importance.keys()
    
    size=len(ProtNames2Importance)
    outmatrix=np.zeros((size,size))
    
    ProtName2Obj={}
    ProtName2GOObj={}
    Prot2Prot_Flow={}
    Prot2Self_Flow={}
    ProtName2Index={}
    ProtIndex=[]
    i=0
    for ProtName in ProtNames2Importance.keys():
        ProtName2Index[ProtName]=i
        ProtIndex.append(ProtName)
        i+=1
        ProtName2Obj[ProtName]=prot_idx['name'][ProtName].single
        ProtName2GOObj[ProtName]=[]
        
        for relationship in ProtName2Obj[ProtName].DESCRIBES:
            ProtName2GOObj[ProtName].append(relationship.start)
    
    for ProtName in ProtName2Obj.keys():
        SelfGOs=ProtName2GOObj[ProtName]
        if not ProtName in Prot2Self_Flow.keys():
            Prot2Self_Flow[ProtName]=0.0
        for GO in SelfGOs:
            #Try to unbuffer this part and see if the performance goes down drastically
            Prot2Self_Flow[ProtName]+=(ProtNames2Importance[ProtName])*GO['inf']
        
        
    for pair in itertools.combinations(ProtName2Obj.keys(),2): 
        # We are not doing self-alignement for each protein present
        CommonGOs=list(set(ProtName2GOObj[pair[0]]) & set(ProtName2GOObj[pair[1]]))
        #Check if the GO object is immutable and can be used as a key in a dictionary 
        if not pair in Prot2Prot_Flow.keys():
            Prot2Prot_Flow[pair]=0.0
        
        for GO in CommonGOs:
            #Try to unbuffer this part and see if the performance goes down drastically
            Prot2Prot_Flow[pair]+=(ProtNames2Importance[pair[0]]+ProtNames2Importance[pair[1]])*GO['inf']
        pairvalue=Prot2Prot_Flow[pair]/(Prot2Self_Flow[pair[0]]+Prot2Self_Flow[pair[1]])
        outmatrix[ProtName2Index[pair[0]]][ProtName2Index[pair[1]]]=pairvalue
        outmatrix[ProtName2Index[pair[1]]][ProtName2Index[pair[0]]]=pairvalue
    
    return [ProtIndex,outmatrix]


def BEA_clustering(numpy_matrix):
    # Assumes a distance matrix 
    # Clusters => performs a BEA algorithm matrix clustering
    # Outputs a clustered matrix and 
    
    col={}
    addflavor=np.diag(np.sum(numpy_matrix, axis=-1))
    for i in range(0,numpy_matrix.shape[0]):
        col[i]=numpy_matrix[:,i]+addflavor[:,i]
     
    itlist=range(0,i+1)
    rand.shuffle(itlist)
    col['zeros']=np.zeros(numpy_matrix[:,0].shape)
    
    finalList=itlist[:2]+['zeros']
    for idx in itlist[2:]:
        EnergyDiff={}
        for indx in range(0,len(finalList)-1):
            id1=indx
            id2=(indx+1)
            EnergyDiff[indx]=np.dot(col[finalList[id1]],col[idx])+np.dot(col[idx],col[finalList[id2]])-np.dot(col[finalList[id1]],col[finalList[id2]])
        insindex=(max(EnergyDiff.iteritems(), key=operator.itemgetter(1))[0]+1)
        finalList.insert(insindex, idx)
    returnmatrix=np.zeros(numpy_matrix.shape)
    
    finalList=finalList[:-1]
    for i in range(0,len(finalList)):
        for j in range(0,len(finalList)):
            returnmatrix[i][j]=numpy_matrix[finalList[i]][finalList[j]]
            
    return [finalList, returnmatrix]
    
def split_afterVal(numpy_matrix):
    # Takes a pre-processed matrix according to BEA and calculates the price of split
    # after a given index
    # Outputs the price of a dict of type {index_before_split:Price_of_split} 
    
    EstimMatrix=numpy_matrix+np.diag(np.sum(numpy_matrix, axis=-1))
    
    SplitValue=[]
    for i in range(0,EstimMatrix.shape[0]-1):
        SplitValue.append(np.dot(EstimMatrix[:,i],EstimMatrix[:,i+1]))
    
    return SplitValue

def find_best_splitters(IndexList,numpy_matrix,success_perc):
    # Transform to dict and output splitting sectors and splitted list
    EstimMatrix=numpy_matrix+np.diag(np.sum(numpy_matrix, axis=-1))
    SplitValue = {}
    for i in range(0,EstimMatrix.shape[0]-1):
        SplitValue[i] = np.dot(EstimMatrix[:,i],EstimMatrix[:,i+1])
    
    Ordered_SV = sorted(SplitValue.iteritems(),key=operator.itemgetter(1))
    BreakSV = []
    for i in range(0,max(int(round(len(Ordered_SV)*0.2)),1)):
        BreakSV.append(Ordered_SV[i][0]+1) # +1 is due to the fact that split breaks before the provided index
    BreakSV.sort()
    arrlist=np.split(np.array(IndexList),BreakSV)
    retlist=[]
    for elt in arrlist:
        retlist.append(list(elt))
    
    return retlist
    
def full_stack_BEA_clustering(ProtName2Imp,filterLvl,verbose=False,normalized=False):
    simObj=[]
    if normalized:
        simObj=rapid_naive_group_pairwise_info_normalized(ProtName2Imp)
    else:
        simObj=rapid_naive_group_pairwise_info(ProtName2Imp) 
    simIdxNames=simObj[0]
    simMatrix=simObj[1]
    clustSimObj=BEA_clustering(simMatrix)
    IdxPermutation=clustSimObj[0]
    clustMat=clustSimObj[1]
    clustIdxNames=[simIdxNames[i] for i in IdxPermutation]
    retour=find_best_splitters(clustIdxNames,clustMat,filterLvl)    
    if verbose:
        print simIdxNames
        print simMatrix
        print clustIdxNames
        print clustMat
        print retour
    
    return retour

def get_chunk_stats(size, defin):
    
    query = select([TableBuilder.UNIPROT_Prot]).order_by(func.random()).limit(size)
    quest1 = session.execute(query)
    rows = quest1.fetchall()
    randomList={}
    
    for row in rows:
        randomList[str(row.uniprotid)]=1.0
    
    statmat = rapid_naive_group_pairwise_info_normalized(randomList)

    Vallist = statmat[1].flatten().tolist()

    statchunks=100.0/float(defin)
    
    statchunklist={}
    
    for line in Vallist:
        print "{0:.2f}".format(line*100)
        statchunkord=round(line*100.0/statchunks)
        if not statchunkord in statchunklist.keys():
            statchunklist[statchunkord]=0
        statchunklist[statchunkord]+=1
    
    for key in statchunklist.keys():
        print key,statchunklist[key]
        
    return statmat

def get_rand_eigenv(size):
       
    query = select([TableBuilder.UNIPROT_Prot]).order_by(func.random()).limit(size)
    quest1 = session.execute(query)
    rows = quest1.fetchall()
    randomList={}
    
    for row in rows:
        randomList[str(row.uniprotid)]=1.0
    
    statmat = rapid_naive_group_pairwise_info_normalized(randomList)
    
    eigvals = np.linalg.eigvalsh(statmat[1])
    
    Vallist = eigvals.flatten().tolist()
    
    return Vallist
    

def get_rand_eigenv_stats(size,iterations):
    Llist=[]
    for i in range(0,iterations):
        Llist+=get_rand_eigenv(50)

    return Llist

#result=get_rand_eigenv_stats(100,200)
#for elt in result:
#    print elt

def load_sec_effect(Name):
    quest = session.query(TableBuilder.DRUG2SECEFF).\
                    join(TableBuilder.MedDRASecEffTerms).\
                    filter(TableBuilder.MedDRASecEffTerms.name==Name).\
                    all()
    
    drugs=[]
    for row in quest:
        drugs.append(row.drugName)
    
    ProtDict={}
    for DrugName in drugs:
        dico=HumanTargets_proc.revMapFile
        if DrugName in dico.keys():
            DrugName=dico[DrugName]
        q = session.query(TableBuilder.DRUG2HUMUNIPROT).\
                    filter(TableBuilder.DRUG2HUMUNIPROT.drugName==DrugName).\
                    all()

        for row in q:
            if '_HUMAN' in row.uniprotid:
                if str(row.uniprotid) not in ProtDict.keys():
                    ProtDict[str(row.uniprotid)]=0
                ProtDict[str(row.uniprotid)]+=1
    
    # Hey, we have to eliminate the random samples of secondary effects!!!!
    
    print len(drugs)
    OrdProtDict=sorted(ProtDict.iteritems(),key=operator.itemgetter(1),reverse=True)
    print len(ProtDict), OrdProtDict
    return ProtDict


def random_Sample(drugs_in_sample):
    
    query = select([TableBuilder.DrugName]).order_by(func.random()).limit(drugs_in_sample)
    quest1 = session.execute(query)
    rows = quest1.fetchall()
    
    drugs=[]
    for row in rows:
        drugs.append(row.name)
        
    ProtDict={}
    for DrugName in drugs:
        dico=HumanTargets_proc.revMapFile
        if DrugName in dico.keys():
            DrugName=dico[DrugName]
        q = session.query(TableBuilder.DRUG2HUMUNIPROT).\
                    filter(TableBuilder.DRUG2HUMUNIPROT.drugName==DrugName).\
                    all()

        for row in q:
            if '_HUMAN' in row.uniprotid:
                if str(row.uniprotid) not in ProtDict.keys():
                    ProtDict[str(row.uniprotid)]=0
                ProtDict[str(row.uniprotid)]+=1
    
    # Hey, we have to eliminate the random samples of secondary effects!!!!
    
#    print 'drugs',len(drugs)
#    OrdProtDict=sorted(ProtDict.iteritems(),key=operator.itemgetter(1),reverse=True)
#    print len(ProtDict), OrdProtDict
    return ProtDict

def RandomDrugSampling(iterations,drugs_in_sample):
    #use numpy here!!!!
    matrix=np.zeros((iterations,1028))
    LableList=[]
    LabelDict={}
    for i in range(0,iterations):
        locdict=random_Sample(drugs_in_sample)
        for key in locdict.keys():
            if key not in LableList:
                LableList.append(key)
                LabelDict[key]=len(LableList)-1
            matrix[i,LabelDict[key]]=locdict[key]
    
    matrix=matrix/drugs_in_sample
    
    mean=np.mean(matrix, axis=0).tolist()
    var=np.var(matrix, axis=0).tolist()
    
    FinDict={}
    for i in range(0,len(LableList)):
        FinDict[LableList[i]]=(mean[i],var[i])
        
    OrdSupdict = sorted(FinDict.iteritems(), key=lambda x:x[1],reverse=True)

    print OrdSupdict 
    return FinDict
 
def filter_randomness_Cheby(ProtDict,drugs_No,filterLevel=0.1):
    trueDict={}
    refdict=UD.Drug_50_100
    for elt in UD.Drug_50_100.keys():
        trueDict[elt]=(refdict[elt][0]*drugs_No,refdict[elt][1]*drugs_No*drugs_No)
    
    NrLikehood={}
    for elt in ProtDict.keys():
        if elt in trueDict.keys():
            NrLikehood[elt] = 1.0-min(trueDict[elt][1]/(ProtDict[elt]-trueDict[elt][0])**2,1.0)
        else:     
            NrLikehood[elt]=1.0
            
    sortedNrLikehood=sorted(NrLikehood.iteritems(),key=operator.itemgetter(1),reverse=True)
    
    outputDict={}
    
    for elt in sortedNrLikehood:
        if elt[1]>1.0-filterLevel:
            print elt[0], "{0:.2f}".format(elt[1]*100)+'%', ProtDict[elt[0]], trueDict[elt[0]][0], trueDict[elt[0]][1]
            outputDict[elt[0]]=elt[1]
    
    return outputDict
        
    
    
    

RealDico=load_sec_effect('Pancreatitis')
Imp_filtered=filter_randomness_Cheby(RealDico,48,0.07)

full_stack_BEA_clustering(Imp_filtered,0.5,verbose=True,normalized=True)

testdict={'RHOC_HUMAN':1.0, 'NGF_HUMAN':1.0, 'NTRK1_HUMAN':1.0}

simple_rapid_naive_group_info(testdict)

#print RandomDrugSampling(100,50) 
 
#get_chunk_stats(500,1000)

#TODO: add importance management to the normalization procedure


##dico={'URE3_MYCTU':1.0,'URE1_MYCTU':1.0,'URE2_MYCTU':1.0,'TPIS_MYCTU':1.0,'TPIS_HUMAN':1.0}
#dico=Targets_Stats.get_eff_drugs_Targets()
#print dico
#start=time()
#full_stack_BEA_clustering(dico,0.2,verbose=True,normalized=True)
#print 'naive fetch and clustering done in:', "{0:.2f}".format(time()-start)
#
## TODO: try an eigenvalue-based approach on splitting, in order to eliminate
## the noize from the random matrix graph compound

neo.shutdown()
print "neo4j is down now"     