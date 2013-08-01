'''
Created on Jul 16, 2013

@author: andrei
'''

from neo4j_Declarations.Graph_Declarator import DatabaseGraph
import copy
from time import time
import pickle
import operator
from random import shuffle
import math
from scipy.stats.kde import  gaussian_kde
import numpy as np
from pylab import plot, hist, show

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

def import_UniprotDict():
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots=pickle.load(pickleDump2)
    pickleDump2.close()
    UniprotDict={}
    for elt in Uniprots:
        node=DatabaseGraph.UNIPORT.get(elt)
        altID=node.ID
        UniprotDict[altID]=(elt,ID2displayName[elt])
        UniprotDict[elt]=altID
    return UniprotDict

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
        else:
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
            vs=acceleratedInsert(GO_structure, accelerationDict, elt)
            visited=visited+vs
        visited=list(set(visited))
        for elt in visited:
            if elt not in TimesReached.keys():
                TimesReached[elt]=0
            TimesReached[elt]+=1
    Fle=file('GO_Informativities.dump','w')
    pickle.dump(TimesReached,Fle)
    Fle2=file('accDict.dump','w')
    pickle.dump(accelerationDict,Fle2)

def analyze_GO_Informativities():
    GO_Infos = pickle.load(file('GO_Informativities.dump','r'))
    srted = sorted(GO_Infos.iteritems(), key=operator.itemgetter(1), reverse=True)
    GO_2_Names=pickle.load(file('GO_names.dump','r'))
    i=0
    for key, val in srted[:500]:
        i+=1
        print i, key, val, GO_2_Names[key]
    return GO_Infos

def load_GO_Informativities():
    GO_Infos = pickle.load(file('GO_Informativities.dump','r'))
    return GO_Infos

def load_GO_Structure():
    GO_structure = pickle.load(file('GO_structure.dump','r'))
    return GO_structure

def load_GO_Accesses():
    GO_Accesses=pickle.load(file('GO.dump','r'))
    return GO_Accesses

def load_accDict():
    accDict=pickle.load(file('accDict.dump','r'))
    return accDict

def acceleratedInsert(GO_structure,accelerationDict, name):
    if name in accelerationDict.keys():
        return accelerationDict[name]
    else:
        cumset=set()
        cumset.add(name)
        for subelt in GO_structure[name][0]:
            cumset.update(acceleratedInsert(GO_structure,accelerationDict, subelt))
        accelerationDict[name]=list(cumset)
        return accelerationDict[name]

def unzip(Supset, GO_2_Names):
    nameList=[]
    for elt in Supset:
        nameList.append(GO_2_Names[elt])
    return nameList

def access_Reactional_Pathways(Uniprot2Vals_Set):
    # Get the molecules from the reactome attached to the uniprot set
    SP_ID2Nodes={}
    for SP_ID in Uniprot2Vals_Set:
        nodeGen=DatabaseGraph.UNIPORT.index.lookup(ID=SP_ID)
        if nodeGen!=None:
            for node in nodeGen: # there can only be one
                SP_ID2Nodes[SP_ID]=node
    SP_ID2ReactNodes={}
    for SP_ID, node in SP_ID2Nodes.iteritems():
        nodeGen=node.bothV("is_same")
        if nodeGen!=None:
            SP_ID2ReactNodes[SP_ID]=[]
            for subnode in nodeGen:
                SP_ID2ReactNodes[SP_ID].append(subnode)
    SP_ID2React_ID={}
    React_ID2React_Names={}
    for SP_ID, nodeList in SP_ID2ReactNodes:
        for node in nodeList:
            gen=node.bothV("is_part_of_pathway")
            if gen!=None:
                for subnode in gen:
                    if subnode.element_type=="Pathway":
                        ID=str(subnode).split('/')[-1][:-1]
                        name=subnode.displayName
                        React_ID2React_Names[ID]=name
                        SP_ID2React_ID[SP_ID]=ID
    return SP_ID2React_ID, React_ID2React_Names
        
    # Get the reactional bindings incoming on those molecules or molecules they are directly related to ()

def convert_SP_to_IDs(SP_List):
    Res_Dict={}
    for name in SP_List:
        generator=DatabaseGraph.UNIPORT.index.lookup(ID=name)
        if generator!=None:
            Res_Dict[name]=[]
            for elt in generator:
                ID=str(elt).split('/')[-1][:-1]
                Res_Dict[name].append(ID)
            if len(Res_Dict[name])>1:
                print 'Error: several references!', name, Res_Dict[name]
            else: 
                Res_Dict[name]=Res_Dict[name][0]          
    return Res_Dict

def specialRatio(Number1,Number2,epsilon=1e-7):
    if abs(Number1)<epsilon and abs(Number2)>epsilon:
        return 0.0
    if abs(Number2)<epsilon and abs(Number1)>epsilon:
        return 100
    if abs(Number2)<epsilon and abs(Number2)<epsilon:
        return 'n.a'
    else:
        return abs(float(Number1)/float(Number1+Number2))*100

def get_GO_Term_occurences(Importance_Dict,flat):
    NamesDict=pickle.load(file('GO_names.dump','r'))
    GO_access=pickle.load(file('GO.dump','r'))[0]
    GO_structure=pickle.load(file('GO_structure.dump','r'))
    GO_Infos = pickle.load(file('GO_Informativities.dump','r'))
    accelerationDict={}
    Associated_GOs={}
    ReverseDict={}
    UP_Dict=import_UniprotDict()
    for key in Importance_Dict.keys():
        toVisit=[]
        toVisit=copy.copy(GO_access[UP_Dict[key][0]])
        visited=[]
        while toVisit!=[]:
            elt=toVisit.pop()
            vs=acceleratedInsert(GO_structure, accelerationDict, elt)
            visited=visited+vs
            visited=list(set(visited))
        Associated_GOs[key]=copy.copy(visited)
        for elt in visited:
            if elt not in ReverseDict.keys():
                ReverseDict[elt]=[]
            ReverseDict[elt].append(key)
    Counter={}
    for UP in Importance_Dict.keys():
        GO_List=Associated_GOs[UP]
        for GO in GO_List:
            if GO not in Counter.keys():
                Counter[GO]=0
            if flat:
                Counter[GO]+=1
            else:
                Counter[GO]+=Importance_Dict[UP]
    Definitive={}
    Definitive_full={}
    cummulLocal=len(Associated_GOs.keys())
    cummulTot=len(GO_access.keys())
    for key, val in Counter.iteritems():
        if val>2:
            exp=float(cummulLocal)*float(GO_Infos[key])/float(cummulTot)
            rat=float(val)/exp
            Definitive[key]=rat
            Definitive_full[key]=(rat, val, exp, NamesDict[key], ReverseDict[key])
            #TODO: Confidence estimation?
    srtd=sorted(Definitive.iteritems(), key=operator.itemgetter(1),reverse=True)
    pdf, log_pdf=get_Tirage_stats()
    pdf2, log_pdf2=get_Dictionnary_Stats(Definitive)
    print  'occurence to expected occurence ratio', '\t', 'occurences','\t', 'expected occurences','\t|\t',
    print  'kernel PDF in sample','\t', 'log-kernel PDF in sample','\t|\t',
    print  'kernel PDF in random sets','\t', 'log-kernel PDF in random sets','\t|\t',
    print  'sample PDF/random set PDF', '\t', 'sample log-PDF/random set log-PDF',
    print  '\t|\t', 'GO Term', '\t', 'List of Targets'
    for key, val in srtd:
        p1=float(pdf2(Definitive_full[key][0]))*100
        p2=float(log_pdf2(math.log(Definitive_full[key][0],10))*100)
        p3=float(pdf(Definitive_full[key][0]))*100
        p4=float(log_pdf(math.log(Definitive_full[key][0],10))*100)
        print  "{0:.2f}".format(Definitive_full[key][0]*100)+'%','\t', Definitive_full[key][1],'\t', "{0:.2f}".format(Definitive_full[key][2]),'\t|\t',
        print  "{0:.2f}".format(p1),'%\t', "{0:.2f}".format(p2),'%\t|\t',
        print  "{0:.2f}".format(p3),'%\t', "{0:.2f}".format(p4),'%\t|\t',
        print  "{0:.2f}".format(specialRatio(p1,p3)), '%\t', "{0:.2f}".format(specialRatio(p2,p4)),
        print  '%\t|\t', Definitive_full[key][3], '\t', Definitive_full[key][4]
    return Associated_GOs, Definitive, Definitive_full
        

def align_names2SP():
    from Utils.UNIPROT_Parser import names_Dict
    from configs import Targets_dict2, Targets_File2
    Fle=file(Targets_File2,'r')
    FileDict={}
    i=0
    print len(names_Dict)
    while True:
        i+=1
        line=Fle.readline()
        if not line:
            break
        if i>3:
            words=line.strip('\n').split('\t')
            FileDict[words[0]]=(1,1,1)
            # FileDict[words[0]]=(float(words[1]),float(words[2].strip()),float(words[3].strip()))
    print len(FileDict)
    Name2SP={}
    for elt in FileDict.keys():
        tp=Targets_dict2[elt]
        if type(tp)==list:
            Name2SP[elt]=tp
        else:
            if len(tp)>1:
                Name2SP[elt]=[names_Dict[tp]]
    return Name2SP

def TouchedIDs():
    Names2SP=align_names2SP()
    valuelist=[]
    for elt in Names2SP.values():
        valuelist=valuelist+elt
    SP_to_IDs=convert_SP_to_IDs(valuelist)
    IDList=[]
    errcount=0
    for valLists in Names2SP.values():
        for val in valLists:
            if val in SP_to_IDs.keys():
                IDList.append(SP_to_IDs[val])
            else:
                errcount+=1
    print '444', len(IDList),errcount,len(valuelist)
    pickle.dump(IDList, file('IDList.dump','w'))
    return IDList
    
#     Uniprot_Dict=import_UniprotDict()
#     i=0
#     j=0
#     final_Dict={}
#     for key, vallist in Name2SP.iteritems():
#         for val in vallist:
#             i+=1
#             if val in Uniprot_Dict.keys():
#                 j+=1
#                 final_Dict[val]=FileDict[key]
#     secDict={}
#     for key, val in final_Dict.iteritems():
#         secDict[key]=-val[2]
#         
#     return final_Dict, secDict

def Tirage(sampleSize, flat, iterations):
    '''
    Pulls at random several proteins of a determined size to estimate 
    the error margins
    '''
    GO_Accesses=pickle.load(file('GO.dump','r'))[0]
    UP_Dict=import_UniprotDict()
    SP_List=list(GO_Accesses.keys())
    RestList=[]
    for i in range(0, iterations):
        shuffle(SP_List)
        ImpDict={}
        for item in SP_List[:sampleSize]:
            ImpDict[UP_Dict[item]]=1.0
        RestList.append(get_GO_Term_occurences(ImpDict,flat)[1:2])
        print '\n<===============================>\n'
    pickle.dump(RestList,file('CompressedStats.dump','w'))
    
def get_Tirage_stats():
    RestList=pickle.load(file('CompressedStats.dump','r'))
    dumpFile=file('dumpFile.csv', 'w')
    ValList=[]
    LogValList=[]
    for elt in RestList:
        for key, val in elt[0].iteritems():
            dumpFile.write(str(val)+'\t'+str(math.log(val,10))+'\n')
            ValList.append(val)
            LogValList.append(math.log(val,10))
    pdf=gaussian_kde(np.asarray(ValList))
    log_pdf=gaussian_kde(np.asarray(LogValList))
    return pdf, log_pdf

def get_Dictionnary_Stats(Dictionary):
    ValList=[]
    LogValList=[]
    for key,val in Dictionary.iteritems():
        ValList.append(val)
        LogValList.append(math.log(val,10))
    pdf=gaussian_kde(np.asarray(ValList))
    log_pdf=gaussian_kde(np.asarray(LogValList))
    return pdf, log_pdf

def get_Uniprot_Subset(List_of_GOs):
    raise NotImplementedError
    
    
# TODO: add the modules for matrix operations over the GO annotation
# TODO: add the propagation of the informativity along different GO Terms
# filtr=['biological_process']
# RelDict, SeedSet=get_GO_access(filtr)
# print 1, time()-init
# get_GO_structure(filtr,SeedSet)
# print 2, time()-init
# get_GO_Informativities()
# align_names2SP()
# FD,SD=align_names2SP()
# get_GO_Term_occurences(SD,True)
# Tirage(48,True,100)
# get_Tirage_stats()
# TouchedIDs()
# => Only about 100 uniprots out of 4000 do not point towards the 