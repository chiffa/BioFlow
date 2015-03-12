'''
Created on Jul 10, 2013

@author: andrei
'''

from PolyPharma.configs2 import Hint_csv
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph

def get_Prot2ProtRels():
    docu=open(Hint_csv,"r")
    LocalRelations={}
    i=0
    while True:
        i+=1
        line=docu.readline()
        if not line:
            break
        if i>1:
            fields=line.split('\t')
            if fields[2]!=fields[3]:
                if fields[2] not in LocalRelations.keys():
                    LocalRelations[fields[2]]=[]
                if fields[3] not in LocalRelations.keys():
                    LocalRelations[fields[3]]=[]
                LocalRelations[fields[3]].append(fields[2])
                LocalRelations[fields[2]].append(fields[3])
    return LocalRelations

def get_Uniprots():
    UPDict={}
    for elt in DatabaseGraph.UNIPORT.get_all():
        ID=str(elt).split('/')[-1][:-1]
        primary=DatabaseGraph.UNIPORT.get(ID)
        UPDict[str(primary.ID).split('_')[0]]=primary
    return UPDict

def cross_ref_HiNT(flush):
    RelationDict=get_Prot2ProtRels()
    UniProtRefDict=get_Uniprots()
    Treated=set()
    i=0
    for key in UniProtRefDict.keys():
        if key in RelationDict.keys():
            Treated.add(key)
            for subkey in RelationDict[key]:
                if subkey in UniProtRefDict.keys() and subkey not in Treated:
                    i+=1
                    print key, subkey
                    if flush:
                        DatabaseGraph.is_interacting.create(UniProtRefDict[key], UniProtRefDict[subkey])
    print i, len(Treated)

