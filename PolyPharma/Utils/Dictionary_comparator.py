'''
Created on Jul 22, 2013

@author: andrei

A set of tool allowing loose matching well-suited for biological applications
'''
import Levenshtein as lv
import itertools
import re
from time import time
import operator

def straighten(line):
    wordlist=re.split(' |,|;|-|\.|/',line)
    wordlist=filter(None, wordlist)
    resstring=''
    for wordset in itertools.combinations(wordlist,len(wordset)):
        for word in wordset:
            resstring=resstring+' '+word.strip()
    return resstring

def segment(line):
    stringList=[]
    wordlist=re.split(' |,|;|-|\.|/',line)
    wordlist=filter(None, wordlist)
    for n in range(min(2,len(wordlist)),len(wordlist)):
        for wordset in itertools.combinations(wordlist,n):
            resstring=''
            for word in wordset:
                resstring=resstring+' '+word.strip()
            resstring.strip()
            stringList.append(resstring)
    return stringList
    
def Get_Unique_Key_Derivatives(Dictionnary):
    UniqueSet=set()
    NonUniqueSet=set()
    CandidateDict={}
    for key in Dictionnary.keys():
        derivatives=segment(key)
        CandidateDict[key]=derivatives
        for elt in derivatives:
            if elt in UniqueSet:
                UniqueSet.remove(elt)
                NonUniqueSet.add(elt)
            else :
                if elt not in NonUniqueSet:
                    UniqueSet.add(elt)
    ExpandedDict={}
    for elt in CandidateDict.keys():
        for subelt in CandidateDict[elt]:
            if subelt in UniqueSet:
                ExpandedDict[subelt]=elt
    return ExpandedDict

def Get_Loose_Term_Matching(keySetDict,StrictDict,Fzmatch):
    init1=time()
    init=time()
    FDict={}
    for key in keySetDict.keys():
        FDict[key.lower().strip()]=key
    SDict={}
    for key in StrictDict.keys():
        SDict[key.lower().strip()]=key
    Answer1={}
    Answer2={}
    Answer3={}
    print 'loaded in:', time()-init
    init=time()
    Remainder=[]
    # First: exact simple Matching
    for key in FDict.keys():
        if key in SDict.keys():
            Answer1[FDict[key]]=SDict[key]
        else:
            Remainder.append(key)
    print 'exact matching in:', time()-init, 'reduced from', len(FDict.keys()),'to',len(Remainder)
    init=time()
    # Second exact pa1rtial matching
    # Symmetric expansion is not done to avoid complexity explosion
    ExpansionDict=Get_Unique_Key_Derivatives(SDict)
    print 'expanded in:', time()-init, 'size went from:', len(SDict), 'to', len(ExpansionDict)
    Remainder2=[]
    init=time()
    for elt in Remainder:
        if elt in ExpansionDict.keys():
            Answer2[FDict[elt]]=SDict[ExpansionDict[elt]]
        else:
            Remainder2.append(elt)
    print 'exact partial matching in:', time()-init, 'reduced from', len(Remainder), 'to', len(Remainder2)
    init=time()
    # Third: loose partial matching
    Remainder3=[]
    for elt in Remainder2:
        Answer3[FDict[elt]]=[]
        for secelt in ExpansionDict.keys():
            if lv.ratio(elt,secelt)>Fzmatch:
                Answer3[FDict[elt]].append((SDict[ExpansionDict[secelt]], elt, secelt, lv.ratio(elt,secelt)))
                Remainder3.append(elt)
        Answer3[FDict[elt]].sort(key=operator.itemgetter(3),reverse=True)
        
    Remainder3=list(set(Remainder2)-set(Remainder2))
    print 'fuzzy partial matching in:', time()-init, 'reduced from', len(Remainder2), 'to', len(Remainder3)
    print time()-init1
    return Answer1, Answer2, Answer3, Remainder3

def align_names2SP():
    from Utils.UNIPROT_Parser import names_Dict
    Fle=file('/home/andrei/workspaces/UCSD/NeflanavirSource.csv','r')
    FileDict={}
    i=0
    print len(names_Dict)
    while True:
        i+=1
        line=Fle.readline()
        if not line:
            break
        if i>3:
            words=line.split('\t')
            FileDict[words[0]]=(words[1],words[2],words[3])
    Ans1, Ans2, Ans3, Rem3 = Get_Loose_Term_Matching(FileDict,names_Dict,0.85)
    for key, val in Ans1.items():
        print 1, key, '!|!', val
    for key, val in Ans2.items():
        print 2, key, '!|!', val
    for key, val in Ans3.items():
        print 3, key
        for elt in val[:10]:
            print '\t', elt[0], '\t|\t', elt[1], '\t|\t', elt[2], '\t|\t', elt[3]
    print Rem3


align_names2SP()