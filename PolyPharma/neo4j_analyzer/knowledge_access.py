'''
Created on Jul 16, 2013

@author: andrei
'''

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

def compute_UniprotDict():
    pickleDump2=file('pickleDump2.dump','r')
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots=pickle.load(pickleDump2)
    pickleDump2.close()
    UniprotDict={}
    for elt in Uniprots:
        node=DatabaseGraph.UNIPORT.get(elt)
        altID=node.ID
        UniprotDict[altID]=(elt,ID2displayName[elt])
        UniprotDict[elt]=altID
    Fle=file('Uniprot_Dict.dump','w')
    pickle.dump(UniprotDict,Fle)
    return UniprotDict

def import_UniprotDict():
    UniprotDict=pickle.load(file('Uniprot_Dict.dump','r'))
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
    GO_IDs={}
    rev_GO_IDs={}
    while seedList!=[]:
        ID=seedList.pop()
        VisitedSet.add(ID)
        LocUpList=[]
        LocRegList=[]
        GONode=DatabaseGraph.GOTerm.get(ID)
        GO_Names[ID]=str(GONode.displayName)
        GO_IDs[ID]=str(GONode.ID)
        rev_GO_IDs[str(GONode.ID)]=ID
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
    Fle3=file('GO_IDs.dump','w')
    pickle.dump(GO_IDs,Fle3)
    Fle4=file('rev_GO_IDs.dump','w')
    pickle.dump(rev_GO_IDs,Fle4)
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
    Reverse_Dict={}
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
                Reverse_Dict[elt]=[]
                TimesReached[elt]=0
            Reverse_Dict[elt].append(key)
            TimesReached[elt]+=1
    Fle=file('GO_Informativities.dump','w')
    pickle.dump(TimesReached,Fle)
    Fle2=file('accDict.dump','w')
    pickle.dump(accelerationDict,Fle2)
    Fle3=file('Reverse_dict.dump','w')
    pickle.dump(Reverse_Dict,Fle3)

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
    accelerationDict=load_accDict()
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
            # DONE: Confidence estimation?
    srtd=sorted(Definitive.iteritems(), key=operator.itemgetter(1),reverse=True)
    pdf, log_pdf=get_Tirage_stats()
    pdf2, log_pdf2=get_Dictionnary_Stats(Definitive)
    Fle=file('GO_Terms_Aboundance_estimation.csv','w')
    print  'occurence to expected occurence ratio', '\t', 'occurences','\t', 'expected occurences','\t|\t',
    Fle.write('occurence to expected occurence ratio \t occurences \t expected occurences \t ')
    print  'kernel PDF in sample','\t', 'log-kernel PDF in sample','\t|\t',
    Fle.write('kernel PDF in sample \t log-kernel PDF in sample \t')
    print  'kernel PDF in random sets','\t', 'log-kernel PDF in random sets','\t|\t',
    Fle.write('kernel PDF in random sets \t log-kernel PDF in random sets \t ')
    print  'sample PDF/random set PDF', '\t', 'sample log-PDF/random set log-PDF',
    Fle.write('sample PDF/random set PDF \t sample log-PDF/random set log-PDF ')
    print  '\t|\t', 'GO Term', '\t', 'List of Targets'
    Fle.write('\t GO Term \t List of Targets \n')
    for key, val in srtd:
        p1=float(pdf2(Definitive_full[key][0]))*100
        p2=float(log_pdf2(math.log(Definitive_full[key][0],10))*100)
        p3=float(pdf(Definitive_full[key][0]))*100
        p4=float(log_pdf(math.log(Definitive_full[key][0],10))*100)
        print  key, '\t', "{0:.2f}".format(Definitive_full[key][0]*100)+'%','\t', Definitive_full[key][1],'\t', "{0:.2f}".format(Definitive_full[key][2]),'\t|\t',
        Sstring="{0:.4f}".format(Definitive_full[key][0])+'\t'+str(Definitive_full[key][1])+'\t'+"{0:.2f}".format(Definitive_full[key][2])+'\t'
        print  "{0:.2f}".format(p1),'%\t', "{0:.2f}".format(p2),'%\t|\t',
        Sstring=Sstring+"{0:.4f}".format(p1/100.0)+'\t'+"{0:.4f}".format(p2/100.0)+'\t'
        print  "{0:.2f}".format(p3),'%\t', "{0:.2f}".format(p4),'%\t|\t',
        Sstring=Sstring+"{0:.4f}".format(p3/100.0)+'\t'+"{0:.4f}".format(p4/100.0)+'\t'
        print  "{0:.2f}".format(specialRatio(p1,p3)), '%\t', "{0:.2f}".format(specialRatio(p2,p4)),
        Sstring=Sstring+"{0:.4f}".format(specialRatio(p1,p3)/100.0)+'\t'+"{0:.4f}".format(specialRatio(p2,p4)/100.0)
        print  '%\t|\t', Definitive_full[key][3], '\t', Definitive_full[key][4]
        Sstring=Sstring+'\t'+str(Definitive_full[key][3])+'\t'+str(Definitive_full[key][4])+'\n'
        Fle.write(Sstring)
    Fle.close()
    return Associated_GOs, Definitive, Definitive_full
        

def align_names2SP():
    from PolyPharma.Utils.UNIPROT_Parser import names_Dict
    # If overington, switch comments on the two next lines
    # from configs import Targets_dict2 as Targets_dict
    # from congigs import Targets_File2 as Targets_File
    from PolyPharma.configs import Targets_dict, Targets_File
    Fle=file(Targets_File,'r')
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
            # FileDict[words[0]]=(1,1,1)
            FileDict[words[0]]=(float(words[1]),float(words[2].strip()),float(words[3].strip()))
    print len(FileDict)
    Name2SP={}
    for elt in FileDict.keys():
        tp=Targets_dict[elt]
        if type(tp)==list:
            Name2SP[elt]=tp
        else:
            if len(tp)>1:
                Name2SP[elt]=[names_Dict[tp]]
    # Comment Out if used with Overinton's data
    Uniprot_Dict=import_UniprotDict()
    i=0
    j=0
    final_Dict={}
    for key, vallist in Name2SP.iteritems():
        for val in vallist:
            i+=1
            if val in Uniprot_Dict.keys():
                j+=1
                final_Dict[val]=FileDict[key]
    secDict={}
    for key, val in final_Dict.iteritems():
        secDict[key]=-val[2]
    return final_Dict, secDict
#     return Name2SP

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

def rebuild():
    filtr=['biological_process']
    RelDict, SeedSet=get_GO_access(filtr)
    get_GO_structure(filtr,SeedSet)
    get_GO_Informativities()
    compute_UniprotDict()

def get_Reference_Flow(GO_Node_ID,GO_ID,py_mongo_collection, GOs2UP_Node_IDs):
    if py_mongo_collection.find({'GO_ID': GO_ID}).count()>0:
        for pointer in py_mongo_collection.find({'GO_ID': GO_ID}):
            Ref_Flow=pickle.loads(pointer['Info_Array'])
            return Ref_Flow
    UPs=GOs2UP_Node_IDs[GO_Node_ID]
    serUPs=pickle.dumps(set(UPs))
    if py_mongo_collection.find({'UP_Set':serUPs}).count()>0:
        for pointer in  py_mongo_collection.find({'UP_Set':serUPs}):
            Ref_Flow=pickle.loads(pointer['Info_Array'])
            post={'GO_ID':GO_ID,'UP_Set':serUPs,'Info_Array':pickle.dumps(Ref_Flow)}
            py_mongo_collection.insert(post)
            return Ref_Flow
    from Matrix_Puller import Compute_and_Store_circulation
    Ref_Flow=Compute_and_Store_circulation(GO_ID,UPs,py_mongo_collection)
    return Ref_Flow

def get_Local_Flow(UP_List,py_mongo_collection):
    serUPs=pickle.dumps(set(UP_List))
    if py_mongo_collection.find({'UP_Set':serUPs}).count()>0:
        for pointer in  py_mongo_collection.find({'UP_Set':serUPs}):
            local_Flow=pickle.loads(pointer['Info_Array'])
            return local_Flow
    from Matrix_Puller import Compute_and_Store_circulation
    local_Flow=Compute_and_Store_circulation('',UP_List,py_mongo_collection)
    return local_Flow

def compare_Local_Info_to_Ref(Go,UP_List):
    '''
    GO is the Node # of the associated GO term
    
    UP_List is a list of Swissprot_IDs reached by a GO term in a particular configuration
    '''
    from PolyPharma.configs import ref_coll,data_coll
    Go2Node_IDs=pickle.load(file('GO_IDs.dump','r'))
    GOs2UP_Node_IDs=pickle.load(file('Reverse_dict.dump','r'))
    SP_ID2Node_ID=import_UniprotDict()
    re_UP_List=[]
    for elt in UP_List:
        re_UP_List.append(SP_ID2Node_ID[elt])
    refFlow=get_Reference_Flow(Go, Go2Node_IDs[Go], ref_coll, GOs2UP_Node_IDs)
    locFlow=get_Local_Flow(re_UP_List,data_coll)
    print refFlow
    print locFlow

def Compute_GO_sp_InfoCirc(Under_N,Over_N):
    GO_UnderN=[]
    GOs2UP_Node_IDs=pickle.load(file('Reverse_dict.dump','r'))
    Go2Node_IDs=pickle.load(file('GO_IDs.dump','r'))
    for GO, val in GOs2UP_Node_IDs.iteritems():
        if len(val)<Under_N and len(val)>Over_N:
            GO_UnderN.append(GO)
    from PolyPharma.configs import ref_coll
    i=0
    l=len(GO_UnderN)
    for GO in GO_UnderN:
        i+=1
        print 'processing :', i, 'out of', l 
        get_Reference_Flow(GO, Go2Node_IDs[GO], ref_coll, GOs2UP_Node_IDs)

def broadcaset_Uniports(Uniprots,UPNode_IDs_2Proteins_IDs_List,NodeID2MatrixNumber,MeanInfos):
    for UP_ID in UPNode_IDs_2Proteins_IDs_List.keys():
        # @attention: patch over here too.
        # TODO: correct the matrix retrieval implementation
        if len(UPNode_IDs_2Proteins_IDs_List[UP_ID])>1:
            print UP_ID, 'problem!!!!'
        else:
            Prot_ID=UPNode_IDs_2Proteins_IDs_List[UP_ID][0]
            if Prot_ID in NodeID2MatrixNumber.keys(): 
                MeanInfos[NodeID2MatrixNumber[UP_ID],0]=MeanInfos[NodeID2MatrixNumber[Prot_ID],0]
    return MeanInfos

def get_Max_Informativities(FilteringAbsolute, FilteringFraction):
    '''
    For each term, computes the set of GO Terms for which a Uniprot term have a maximal informativity
    '''
    from PolyPharma.configs import ref_coll
    from Matrix_Puller import load_Uniprot_Attachments
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    UPNode_IDs_2Proteins_IDs_List=load_Uniprot_Attachments()
    print len(NodeID2MatrixNumber.keys()) 
    # Broadcast the uniprot attachments:
    GO2Column={}
    Column2GO={}
    Lst=[]
    Lst2=[]
    GO2TotalFlow={}
    i=0
    for InfArrayDict in ref_coll.find():
        print 'trating:', InfArrayDict['GO_ID'], i
        TrueArray = broadcaset_Uniports(Uniprots,UPNode_IDs_2Proteins_IDs_List,NodeID2MatrixNumber,pickle.loads(InfArrayDict['Info_Array']))
        TrueSet = (InfArrayDict['GO_ID'],pickle.loads(InfArrayDict['UP_Set']),TrueArray)
        GO2Column[TrueSet[0]] = i
        Column2GO[i] = TrueSet[0]
        GO2TotalFlow[TrueSet[0]] = len(TrueSet[1])
        Lst.append(TrueSet[2]/len(TrueSet[1]))
        Lst2.append(TrueSet[2]/float(len(TrueSet[1])**2/2.0))
        i+=1
    # Use term informativity to compute the basis of flow
    Mat1=np.concatenate(tuple(Lst),axis=1)
    Mat2=np.concatenate(tuple(Lst2),axis=1)
    print Mat1.shape
    print Mat2.shape
    # For each Uniprot, get the highest absolute Flow
    # and highest percentage flow
    HighestAbsolute_Dict={}
    HighestRelative_Dict={}
    
    for UP_ID in Uniprots:
        row_Num=NodeID2MatrixNumber[UP_ID]
        besthit1=np.argmax(Mat1[row_Num,:])
        besthit1_val=Mat1[row_Num,besthit1]
        besthit2=np.argmax(Mat2[row_Num,:])
        besthit2_val=Mat2[row_Num,besthit2]

        if besthit1_val>FilteringAbsolute:
            HighestAbsolute_Dict[UP_ID]=(Column2GO[besthit2],besthit1_val)
        if besthit2_val>FilteringFraction:
            HighestRelative_Dict[UP_ID]=(Column2GO[besthit1],besthit2_val)
    
    print len(Uniprots), len(HighestAbsolute_Dict.keys()), len(HighestRelative_Dict.keys())
    
    NamesDict=pickle.load(file('GO_names.dump','r'))
    Rev_GO_IDs=pickle.load(file('rev_GO_IDs.dump','r'))
    
    for key,val in HighestRelative_Dict.iteritems():
        print ID2displayName[key],'\t', val[0],'\t', NamesDict[Rev_GO_IDs[val[0]]],'\t', val[1]
    
    print '\n\n\n<===========================================================>\n\n\n'
    
    for key,val in HighestAbsolute_Dict.iteritems():
        print ID2displayName[key],'\t', val[0],'\t', NamesDict[Rev_GO_IDs[val[0]]],'\t', val[1]
    
    # Reincode as arrays
    
    Max_GOs = np.chararray((len(MatrixNumber2NodeID.keys()),1),itemsize=200)
    Max_GOs[:]=''
    Max_GO_Names = np.zeros((len(MatrixNumber2NodeID.keys()),1))
    Max_GO_val = np.zeros((len(MatrixNumber2NodeID.keys()),1))
    
    for key,val in HighestRelative_Dict.iteritems():
        index=NodeID2MatrixNumber[key]
        Max_GOs[index,0]=str(NamesDict[Rev_GO_IDs[val[0]]])
        Max_GO_Names[index,0]=val[0]
        Max_GO_val[index,0]=val[1]
    
    for key,val in HighestAbsolute_Dict.iteritems():
        index=NodeID2MatrixNumber[key]
        Max_GOs[index,0]=str(NamesDict[Rev_GO_IDs[val[0]]])+' *'
        Max_GO_Names[index,0]=val[0]
        Max_GO_val[index,0]=val[1]/float(GO2TotalFlow[val[0]])*2.0
    
    finmatrix=np.concatenate((Max_GOs,Max_GO_Names,Max_GO_val),axis=1)
    Fle=file('finmatrix.dump','w')
    pickle.dump(finmatrix,Fle)
    return HighestAbsolute_Dict, HighestRelative_Dict, finmatrix
        
                
        
    
# rebuild()
# Tirage(48,True,100)
# FD,SD = align_names2SP()
# get_GO_Term_occurences(SD,True)
# get_Tirage_stats()
# TouchedIDs()
# Compute_GO_sp_InfoCirc(10,3)
# get_Max_Informativities(1.5,0.35)







# => Only about 100 uniprots out of 4000 do not point towards the GOs


    # Step 1: 
    
    # if any given uniprot routes over 30% of information related to this flow 
            
    # In fact we are only interested in the General Reference computation. 
    # The only reason we limit ourselfs to N terms is because we cannot afford computation of all the GO terms.
    # => Can we recover all the GO terms that have < 50 affected uniprots in our set and if there are not too much of them, compute the highest impact ones?
    # the use here is for the Overingtonicity.
    
    # Then, we can take all the targets affected by a drug and perform a computation ponderated by the 
    
    
    
    # recover the absolute value of the difference between the two in percentage of the original value
    #    => Euh... ce sera egale a la valeur de nombre des uniprots de difference, non?
    
    # recover the most important terms for the flow (optionally)
    
    
    # check the availability of the existing reference GO, if not of the UP Set associated to it
        # Lookup by GO ID
        # If fails, lookup by Protein association
        
        # If fails, compute and insert
        
    # for the imported GO
        # check the availability of the Prot_Set
        # If fails, compute and insert, by using only Prot_Set as key
    
# DONE: add UniprotLists Accessed for different GOs
    
    
# CANCEL: add the modules for matrix operations over the GO annotation
# DONE: add the propagation of the informativity along different GO Terms

# CANCEL: possible acceleration with by using the reverse_dicitonary filter: count how many targets are present in the filter.
# However the use of a fully accelerated Dict is pretty much as efficient.


# Now, it is pretty interesting to notice that anything even remotedly related to the central neural system is really badly treated by this method.
# Same thing for the 
