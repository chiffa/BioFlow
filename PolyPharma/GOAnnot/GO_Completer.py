'''
Created on Mar 7, 2013

@author: akucahravy
'''

from DBLoader import TableBuilder
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from AcceleratedDB import GO_Rapidfire
#from AcceleratedDB import UNIPROT_FullGO
from DBLoader.TableBuilder import UNIPROT_FullGO2
from time import time
import logging
import math
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='/home/akucahravy/DB/Annot_Density.txt',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

regTerms=['negatively_regulates','positively_regulates','regulates','negatively_regulated_by','positively_regulated_by','regulated_by']
isaTerms=['is_a','part_of']

def ReciprGOTerms():
    '''
    This method ensures that all the GO relations stored in the GO_Relation
    table are symmetrical regarding the regulatory relations
    '''
    allRelations=session.query(TableBuilder.GO_Relation).\
                        filter(TableBuilder.GO_Relation.reltype!='is_a').\
                        filter(TableBuilder.GO_Relation.reltype!='part_of').\
                        all()
    
    for relation in allRelations:
        if relation.reltype=='negatively_regulates':
            newRel=TableBuilder.GO_Relation(relation.goto, 'negatively_regulated_by', relation.gofrom)
            if newRel not in allRelations:
                session.add(newRel)
        if relation.reltype=='positively_regulates':
            newRel=TableBuilder.GO_Relation(relation.goto, 'positively_regulated_by', relation.gofrom)
            if newRel not in allRelations:
                session.add(newRel)
        if relation.reltype=='regulates':
            newRel=TableBuilder.GO_Relation(relation.goto, 'regulated_by', relation.gofrom)
            if newRel not in allRelations:
                session.add(newRel)
    session.commit()

         
    
def recoverRelations(selectedGOTerm):
    '''
    Recovers the relation for a selected GO term 
    '''
    allFromTerms=session.query(TableBuilder.GO_Relation).\
                            filter(TableBuilder.GO_Relation.gofrom==selectedGOTerm).\
                            all()
    mapp={}
    for elemnt in allFromTerms:
        mapp[str(elemnt.goto)]=str(elemnt.reltype)
    
    return mapp

def getGORapid(GOTerm,timeout=3):
    '''
    Implements a recursive buffering to recover the GO annotation for a given protein
    Attention, this method is however very slow and have to be run only on the minimal
    set of proteins, for it takes a large time for completion
    
    Attention, this method is in-depth first and is redundant. It might be improved by in-breadth
    first approach
    '''
    
    logging.debug('%s >>>>>>>>>>>>>>>>>>>>>>>>\t\t\t %s \t %s', time(), GOTerm, timeout)
    
    subtime = time()
    
    if session.query(GO_Rapidfire).\
                filter(GO_Rapidfire.rootGO==GOTerm).\
                filter(GO_Rapidfire.timeout==timeout).\
                count()>0:
        
        reslist=session.query(GO_Rapidfire).\
                filter(GO_Rapidfire.rootGO==GOTerm).\
                filter(GO_Rapidfire.timeout==timeout).\
                all()
                
        DerivList=[]
        
        for elt in reslist:
            DerivList.append(str(elt.otherGO))
        
        logging.debug( '%s <<<< via DB in:         \t%s\t\t%s\t%s',time(), "{0:.4f}".format(time()-subtime), GOTerm, timeout)

        return DerivList
    
    else:
        GODict=recoverRelations(GOTerm)
        FullRes=[]
        
        for GOTrm in GODict.keys():
            if GODict[GOTrm] in regTerms:
                if timeout>0:
                    partRes=getGORapid(GOTrm,timeout-1)
                    FullRes=FullRes+[GOTrm]+partRes
                    
            else:
                partRes=getGORapid(GOTrm,timeout)
                FullRes=FullRes+[GOTrm]+partRes
                
        FullRes=list(set(FullRes))
            
        for GOTiterm in FullRes:
            session.add(GO_Rapidfire(GOTerm, timeout, GOTiterm))
            session.commit()
        
        logging.debug('%s <<<< via exploration in:\t%s\t\t%s\t%s',time(), "{0:.4f}".format(time()-subtime), GOTerm, timeout)
        
        return FullRes
    
    

def exploreUniprotRapid(UNIPROTid,timeout):
    '''
    Builds and recovers the full GO annotation Tree for one single UNIPROT id
    
    attention :: there are two tables being filled synchroniousely, due to a poor
            DB structure design and a need for performing extensive joins
    '''
    
    if session.query(UNIPROT_FullGO2).\
                filter(UNIPROT_FullGO2.UNIPROTID==UNIPROTid).\
                filter(UNIPROT_FullGO2.timeout==timeout).\
                count()>0:
        reslist=session.query(UNIPROT_FullGO2).\
                filter(UNIPROT_FullGO2.UNIPROTID==UNIPROTid).\
                filter(UNIPROT_FullGO2.timeout==timeout).\
                all()
        DerivList=[]
        for elt in reslist:
            DerivList.append(str(elt.GOs))
        logging.info('recovered')
        return DerivList
    
    else:
        allGOTerms=session.query(TableBuilder.UNIPROT2GO).\
                            filter(TableBuilder.UNIPROT2GO.uniprotid==UNIPROTid).\
                            all()
        result=[]

        for row in allGOTerms:
            result.append(str(row.goid))
    
        resultList=result[:]

        for term in result:
            resultList =resultList+getGORapid(term,timeout)
        
        resultList=list(set(resultList))
        
        for termiterm in resultList:
            session.add(UNIPROT_FullGO2(UNIPROTid, timeout, termiterm))
         
        logging.info('done anew')
        return resultList

def CoverDB(genTimeout=3):
    '''
    Covers the whole Uniprot Database by checking if a UNIPROT was already annotated 
    and if there is still non-annotated entries, it annotates them
    
    genTimeout is the total timeout of coverage (how may "regulates" relations
    would be passed before the method accepts no more "regulatory" statements)
    
    An empiric value of 3 seems to be pretty wella adapted.
    
    '''
    Query=session.query(TableBuilder.UNIPROT_Prot).all()
    
    UNIPROTids=[]
    
    for row in Query:
        UNIPROTids.append(str(row.uniprotid))
    
    logging.info('%s UniprotIDS to proceed', len(UNIPROTids))
    
    Klen=len(UNIPROTids)
    Klen=Klen/100
    print Klen
    
    i=0
    
    for iD in UNIPROTids:
        i=i+1
        logging.info('%s %s processing', i, iD)
        start = time()
        logging.info( str(exploreUniprotRapid(iD,genTimeout)))
        if i%Klen==0:
            session.commit()
        elapsed = (time() - start)
        logging.info('>>>> processed in: %s',elapsed)
    
        session.commit()


def InformativityTerms(Timeout=3):
    '''
    Computes the informativity of all the terms and returns it to the 
    TABLE
    '''
    
    GO2Count=[]
    
    Query=session.query(TableBuilder.GO_Term).all()
    
    for row in Query:
        GO2Count.append(str(row.termid))
    
    RootTerms={'biological_process': math.log(session.query(UNIPROT_FullGO2).filter(UNIPROT_FullGO2.GOs=='0008150').count()+1),
               'cellular_component': math.log(session.query(UNIPROT_FullGO2).filter(UNIPROT_FullGO2.GOs=='0005575').count()+1),
               'molecular_function': math.log(session.query(UNIPROT_FullGO2).filter(UNIPROT_FullGO2.GOs=='0003674').count()+1)
               }
    
    for GoId in GO2Count:
        
        q=session.query(UNIPROT_FullGO2).\
                    filter(UNIPROT_FullGO2.GOs==GoId).\
                    filter(UNIPROT_FullGO2.timeout==Timeout).\
                    count()
        
        
        if q>0:
            GoIdType=session.query(TableBuilder.GO_Term).\
                    filter(TableBuilder.GO_Term.termid==GoId).\
                    filter(UNIPROT_FullGO2.timeout==Timeout).\
                    one().namespace
            
            print GoId, GoIdType ,RootTerms[GoIdType]/math.log(q+1)
            session.add(TableBuilder.GO_Informativity(GoId,RootTerms[GoIdType]/math.log(q+1)))
    
    session.commit()
    
def analyzeAnnotationDensity(Timeout=3):        
    
    Query=session.query(TableBuilder.UNIPROT_Prot).all()
    
    UNIPROTids=[]
    
    for row in Query:
        UNIPROTids.append(str(row.uniprotid))
        
    for UNIPTID in UNIPROTids:
        query1=session.query(TableBuilder.GO_Informativity).\
                        join(TableBuilder.GO_Term).\
                        join(TableBuilder.UNIPROT_FullGO2).\
                        filter(TableBuilder.UNIPROT_FullGO2.UNIPROTID==UNIPTID).\
                        filter(UNIPROT_FullGO2.timeout==Timeout).\
                        all()
        
        #print 'debug: ', UNIPTID, query1
    
        Sum=0.0
        
        for row in query1:
            Sum+=math.log(row.informativity)
            
        logging.info('%s \t %s \t %s', UNIPTID, Sum, len(query1))
        
        session.add(TableBuilder.UNIPROT_SpecInf(UNIPTID,Timeout,Sum,len(query1)))
    
    session.commit()
            
    
#<==================================================>

#CoverDB()

#InformativityTerms()
    
#analyzeAnnotationDensity()  