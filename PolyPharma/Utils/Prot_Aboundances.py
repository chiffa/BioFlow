'''
Created on Jul 26, 2013

@author: andrei
'''

import pickle

# TODO: input the aboundances into the database, relative to the uniprot terms

def compute_Prot_Aboundances():
    '''
    performs a computationally expensive part of protein aboundance computation
    and stores the result in a pickle file
    '''
    from PolyPharma.Utils.UNIPROT_Parser import get_access_dicts
    from PolyPharma.configs import Prot_abound

    access_dict = get_access_dicts()
    SP2Aboundances={}
    Fle=file(Prot_abound,'r')
    i=0
    while True:
        line=Fle.readline()
        if not line:
            break
        print line
        if line[0]!='#':
            words=line.split('\t')[1:3]
            if  words[0].split('.')[1] in access_dict.keys():
                i+=1
                print words,
                print words[0].split('.')[1],
                print access_dict[words[0].split('.')[1]]
                SP_ID=access_dict[words[0].split('.')[1]]
                SP2Aboundances[SP_ID]=float(words[1])
    FullFle=file('SP2Aboundaces.dump','w')
    pickle.dump(SP2Aboundances,FullFle )
    FullFle.close()
    print i
    return SP2Aboundances

def load_Prot_Aboundances():
    '''
    rapidly loads the protein aboundance data stored in the pickle file
    '''
    SP2Aboundances=pickle.load(file('SP2Aboundaces.dump','r'))
    return SP2Aboundances

def load_Prot_Aboundances_NodeIDs():
    SP2Aboundances=pickle.load(file('../Utils/SP2Aboundaces.dump','r'))
    ID2Aboundances={}
    from PolyPharma.neo4j_analyzer.knowledge_access import convert_SP_to_IDs
    SP2IDs=convert_SP_to_IDs(SP2Aboundances.keys())
    for key,val in SP2Aboundances.iteritems():
        if key in SP2IDs.keys():
            ID2Aboundances[SP2IDs[key]]=val
    return ID2Aboundances

# compute_Prot_Aboundances()
# TODO: move the information onto the Uniprot Nodess
ID2Aboundances=load_Prot_Aboundances_NodeIDs()