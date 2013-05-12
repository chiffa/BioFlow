'''
Created on Apr 5, 2013

@author: akucahravy

Generic class build to contain operations on the numpy matrixes related to 
clustring, normalization or eigenvalue/eignevector retrieval
'''
from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import func, select
import Neo4j_Transactions as N4jT 
import operator
import numpy as np
import random as rand

np.set_printoptions(precision=2)

lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb', echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

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
        simObj=N4jT.rapid_naive_group_pairwise_info_normalized(ProtName2Imp)
    else:
        simObj=N4jT.rapid_naive_group_pairwise_info(ProtName2Imp) 
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
    
    statmat = N4jT.rapid_naive_group_pairwise_info_normalized(randomList)

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
    
    statmat = N4jT.rapid_naive_group_pairwise_info_normalized(randomList)
    eigvals = np.linalg.eigvalsh(statmat[1])
    Vallist = eigvals.flatten().tolist()
    
    return Vallist
    

def get_rand_eigenv_stats(size,iterations):
    # Perform batch execution to avoid switching on/off the neo4j database
    # sonstantly
    
    Llist=[]
    for i in range(0,iterations):
        Llist+=get_rand_eigenv(50)

    return Llist

#result=get_rand_eigenv_stats(100,200)
#for elt in result:
#    print elt