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
import configs as conf

def straighten(line):
    wordlist=re.split(' |,|;|-|\.|/|:|(|)',line)
    wordlist=filter(None, wordlist)
    if 'protein' in wordlist:
        wordlist.remove('protein')
    resstring=''
    for wordset in itertools.combinations(wordlist,len(wordlist)):
        for word in wordset:
            resstring=resstring+' '+word.strip()
        resstring=str(resstring.strip())
    return resstring

def segment(line):
    stringList=[]
    wordlist=re.split(' |,|;|-|\.|/|:|(|)',line)
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

def Get_Loose_Term_Matching(keySetDict,StrictDict,Fzmatch1,Fzmatch2):
    # TODO: entropy-based algorithm?
    # TODO: abbreviation management?
    # TODO: permutations within the names?
    # TODO: add greek names conversion
    init1=time()
    init=time()
    FDict={}
    for key in keySetDict.keys():
        FDict[straighten(key.lower().strip())]=key
    SDict={}
    for key in StrictDict.keys():
        SDict[straighten(key.lower().strip())]=key
    Answer1={}
    Answer2={}
    Answer2_5={}
    Answer3={}
    print 'loaded in:', time()-init
    init=time()
    Remainder=[]
    # First: exact simple and exact "proteinized' matching (geneName+protein) 
    for key in FDict.keys():
        if key in SDict.keys():
            Answer1[FDict[key]]=SDict[key]
        else:
            Remainder.append(key)
    print 'exact matching in:', time()-init, 'reduced from', len(FDict.keys()),'to',len(Remainder)
    init=time()
    # Second exact partial matching
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
    # Third: loose full matching 
    Remainder2_5=[]
    for elt in Remainder2:
        Answer2_5[FDict[elt]]=[]
        for secelt in SDict.keys():
            if  lv.ratio(elt,secelt)>Fzmatch1:
                Answer2_5[FDict[elt]].append((SDict[secelt], elt, secelt, lv.ratio(elt,secelt)))
                Remainder2_5.append(elt)
        if len(Answer2_5[FDict[elt]])<1:
            del Answer2_5[FDict[elt]]
        else:
            Answer2_5[FDict[elt]].sort(key=operator.itemgetter(3),reverse=True)
    Remainder2_5=list(set(Remainder2)-set(Remainder2_5))    
    print 'fuzzy full matching in:', time()-init, 'reduced from', len(Remainder2), 'to', len(Remainder2_5)
    # Fourth: loose partial matching
    Remainder3=[]
    for elt in Remainder2_5:
        Answer3[FDict[elt]]=[]
        for secelt in ExpansionDict.keys():
            if lv.ratio(elt,secelt)>Fzmatch2:
                Answer3[FDict[elt]].append((SDict[ExpansionDict[secelt]], elt, secelt, lv.ratio(elt,secelt)))
                Remainder3.append(elt)
        if len(Answer3[FDict[elt]])<1:
            del Answer3[FDict[elt]]
        else:
            Answer3[FDict[elt]].sort(key=operator.itemgetter(3),reverse=True)
    Remainder3=list(set(Remainder2_5)-set(Remainder3))
    Remainder4=[]
    for elt in Remainder3:
        Remainder4.append(FDict[elt])
    print 'fuzzy partial matching in:', time()-init, 'reduced from', len(Remainder2_5), 'to', len(Remainder3)
    print time()-init1
    return Answer1, Answer2, Answer2_5, Answer3, Remainder4

def align_names2SP():
    from Utils.UNIPROT_Parser import names_Dict
    Fle=file(conf.Drug_Targets,'r')
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
            FileDict[words[0]]=('')
    Ans1, Ans2, Ans2_5, Ans3, Rem3 = Get_Loose_Term_Matching(FileDict,names_Dict,0.80,0.80)
    print len(Ans1.keys()), len(Ans2.keys()), len(Ans2_5.keys()), len(Ans3.keys()), len(Rem3)
    
    print '{'
    for key, val in Ans1.items():
        print '\''+key+'\'', ':', '\''+val+'\''
    print '}'
    for key, val in Ans2.items():
        print 2, key, '\t!|!\t', val
    for key, val in Ans2_5.items():
        print 3, key
        for elt in val[:10]:
            print '\t', elt[0], '\t!|!\t', elt[1], '\t!|!\t', elt[2], '\t!|!\t', elt[3]
    for key, val in Ans3.items():
        print 4, key
        for elt in val[:10]:
            print '\t', elt[0], '\t!|!\t', elt[1], '\t!|!\t', elt[2], '\t!|!\t', elt[3]
    print Rem3

align_names2SP()