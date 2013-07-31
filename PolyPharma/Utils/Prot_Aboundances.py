'''
Created on Jul 26, 2013

@author: andrei
'''

import pickle

def compute_Prot_Aboundances():
    '''
    performs a computationally expensive part of protein aboundance computation
    and stores the result in a pickle file
    '''
    from UNIPROT_Parser import access_dict
    from configs import Prot_abound 
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

# compute_Prot_Aboundances()
SP2Aboundances=load_Prot_Aboundances()