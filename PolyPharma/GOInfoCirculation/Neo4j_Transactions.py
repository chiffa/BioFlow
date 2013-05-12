'''
Created on Apr 5, 2013

@author: akucahravy

 This module was build to hold all the classes responsible for the communication
 Neo4j database hosted on the local machine, so that the calls from the otehr classes
 could be brought to the minimum.

'''

import os
os.environ['NEO4J_PYTHON_JVMARGS']='-Xms128M -Xmx512M'
os.environ['JAVA_HOME']='/usr/lib/jvm/java'
from neo4j import GraphDatabase

from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from time import time
import datetime
import math
import operator
import itertools
import numpy as np

np.set_printoptions(precision=2)

import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

 
neo=GraphDatabase('/home/akucahravy/DB/neo4j')



#If you need to create a database first, uncomment the following paragrap
#and comment the one that follows

#with neo.transaction:
#    proteins=neo.node()
#    gos=neo.node()
#    neo.reference_node.PROTEINS(proteins)
#    neo.reference_node.GOS(gos)
#    prot_idx=neo.node.indexes.create('proteins')
#    gos_idx=neo.node.indexes.create('gos')

# Recover a pre-existing database     
with neo.transaction:
    prot_idx=neo.node.indexes.get('proteins')
    gos_idx=neo.node.indexes.get('gos')
    for rel in neo.reference_node.GOS:
        gos = rel.end
    for rel in neo.reference_node.PROTEINS:
        proteins = rel.end

#Caches for proteins, GOs and edges between them for the whole insertion procedure
ProtNamesBuffer=set() # Protein Cache
GONamesBuffer=set() # GO Cache
EdgesBuffer=set() # Cache for the edges linking the two

def create_Protein(name,Importance):
    '''
    Transactional single Protein node creation
    '''
    if prot_idx['name'][name].single!=None:
        return prot_idx['name'][name].single
    else:
        with neo.transaction:
            ProtObj = neo.node(name=name,imp=Importance)
            ProtObj.INSTANCE_OF(proteins)
            prot_idx['name'][name]=ProtObj 
        return ProtObj

def create_GO(GOId, Informativity):
    '''
    Transactional single GO node creation
    '''
    if gos_idx['name'][GOId].single!=None:
        return gos_idx['name'][GOId].single
    else:
        with neo.transaction:
            GOObj = neo.node(name=GOId,inf=Informativity)
            GOObj.INSTANCE_OF(gos)
            gos_idx['name'][GOId]=GOObj 
        return GOObj

def link(ProtName,GOName):
    '''
    Transactional single Protein node creation
    '''
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
    '''
    Performs a transactional batch creation of nodes and edges with
    prevention of duplication
    
    Parameters: Edge List: a list of edges to be created under the form
                            {(Protein_ID,GO_ID):(Protein_Importance,GOName_informativity)}
    '''
    Prots_TBC={} # list of proteins to be created 
    GOs_TBC={} # list of GOs to be created
    Edges_TBC=[] # edges to be created
    Edges_TC=[] # edges to check
    
    for key in EdgeList.keys():
        if key[0] not in ProtNamesBuffer:               # First use the global buffer (faster then neo4j transaction)
            if prot_idx['name'][key[0]].single==None:   # then check for the presence in neo4j database
                Prots_TBC[key[0]]=EdgeList[key][0]      # if both fail, the protein node have to be created
                Edges_TBC.append(key)                   # Add the edge directly to the "to be added liset", for they cannot exist if one of the nodes it is attaches to is inexistent.
            else:
                Edges_TC.append(key)                # if the node is in the neo4j but not in the global buffer, add it 
                ProtNamesBuffer.add(key[0])         # to speed up the process if it is encountered again during the same insertion session
        else:
            Edges_TC.append(key)                # and in case the protein is in the the database (buffer or not) add the edge that contains it to the "to check" list
        if key[1] not in GONamesBuffer:                 
            if gos_idx['name'][key[1]].single==None:    # Same thing as above, but for the GO object
                GOs_TBC[key[1]]=EdgeList[key][1]
                Edges_TBC.append(key)
            else:
                Edges_TC.append(key) 
                GONamesBuffer.add(key[1])
        else:
            Edges_TC.append(key)
            
    Edges_TC=set(Edges_TC)
    Edges_TBC=set(Edges_TBC)
    Edges_TC=Edges_TC-Edges_TBC     # Just for safety :D
    Edges_TC=Edges_TC-EdgesBuffer   # remove from the edges to be checked all the edges that are already in the edge buffer
    Edges_TC=list(Edges_TC)
    Edges_TBC=list(Edges_TBC)
    
    # Control the Edges_TC list to check if the edges it contains are already present in the neo4j database
    # Group according to the ProtObj,i.e. first value to decrease the number of calls to neo4j
    ProtNames=[]                        # First retrieve the protein objects that are encountered in the Edges_TC list
    for entry in Edges_TC:
        ProtNames.append(entry[0])
    
    ProtNames=list(set(ProtNames))
    
    for ProtName in ProtNames:
        sublist=[]                              # a sub-list for all the GOs already associated to a given protein node
        for key in Edges_TC:            
            if key[0]==ProtName:
                sublist.append(key)
        ProtObj=get_ProtObj(ProtName)
        
        GoIds=[]
        for key in sublist:
            GoIds.append(key[1])
        
        for relationship in ProtObj.DESCRIBES:  # fill it by iterating through proteins within neo4j
            nom=relationship.start['name']
            if nom in GoIds:
                EdgesBuffer.add((ProtName,nom)) # if edge is found within neo4j, add it to the session buffer
    
    Edges_TC=list(set(Edges_TC)-EdgesBuffer)        #clean again the Edges_TC list
    Edges_TBC=list(set(Edges_TC)|set(Edges_TBC))    # fuse it with Edges_TBC
    
    #And finally perform the transactional batch insert!!!!
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
    '''
    Recovers a neo4j protein object provided it's id (unique swissprot name)
    '''
    return prot_idx['name'][ProtName].single

def get_GOObj(GOId):
    '''
    Recovers a neo4j GO object provided it's id (unique short name)
    '''
    return gos_idx['name'][GOId].single

def naive_copy_N_entries(number):
    '''
    A naive transactional insertion, with possible duplication of nodes and 
    links
    
    Input: number => number of first rows that will be copied from the database
     
    '''
    # Recover the first "number" rows from the sqlite database
    Query=session.query(TableBuilder.UNIPROT_FullGO2).\
                            limit(number).\
                            all()
    # Copy to a more convinient format                         
    Prots=[]
    GOs=[]
    links=[]
    for row in Query:
        Prots.append(str(row.UNIPROTID))
        GOs.append(str(row.GOs))
        links.append((str(row.UNIPROTID),str(row.GOs))) 
    
    # Remove duplicates
    Prots=list(set(Prots))
    GOs=list(set(GOs))
    links=list(set(links))
    
    #insert one by one into the neo4j database!
    for ProtName in Prots:
        create_Protein(ProtName, 1.0)
    for GOId in GOs:
        create_GO(GOId, 1.0)
    for Llink in links:
        link(Llink[0],Llink[1])        

def copy_from_to_entries(offset,chunk_size):
    '''
    Batch-copies a provided number of links encoded by the sql database,
    with a possibility to jump several chunks before starting the copying process
    
    Arguments:      offset - number of chunks to jump before copying a chunk
                    chunk_size - number of lines within a chunk (optimal is ~1000 )
    '''
    #read from sql daatabase
    Query=session.query(TableBuilder.UNIPROT_FullGO2).\
                            offset(offset*chunk_size).\
                            limit(chunk_size).\
                            all()
    #transform into a form readable by batch insert
    links={}
    for row in Query:
        links[(str(row.UNIPROTID),str(row.GOs))]=(1.0,1.0)
    # use batch_create method to create all of the links
    batchCreate_links(links)

def load_all_chunks(chunk_size):
    '''
    Batch copies all the links and nodes encoded by a sql database, processing with 
    chunks of size chunk_size.
    
    Arguments:    chunk_size - size of chunks that are being copied.
    '''
    # determine number of chunks:
    Query=session.query(TableBuilder.UNIPROT_FullGO2).count()
    number_of_chunks=Query/chunk_size
    print 'total chunks:', number_of_chunks
    estimatedTime=100
    # iterate chunk insertion while increasing the offset!
    for i in range(0,number_of_chunks+1):
        start=time()                        #  we also measure time that it takes to insert all of the links
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

neo.shutdown()
print "neo4j is down now"  