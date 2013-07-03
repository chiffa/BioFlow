'''
Created on Jul 2, 2013

@author: andrei
'''

import configs as conf

GO_Terms={} 
'''
Parse dictionary is of the type GO_Term_ID: {}
'''

def fill_GO_Terms():
    docu=open(conf.GeneOntology,"r")
    ## other_types=set() # keeps the types of the "other" relations
    # ['subset', => ignore
    # 'comment', => ignore
    # 'exact_synonym',
    # 'consider', => ignore
    # 'relationship', => integrate
    #    'part_of'
    #    'regulates'
    #    'positively_regulates'
    #    'negatively_regulates'
    # 'related_synonym',
    # 'narrow_synonym',
    # 'broad_synonym',
    # 'replaced_by',
    # 'alt_id',
    # 'xref_analog']
    localDictionary={}
    i=0
    blocks=0
    Block=False
    Obsolete=False
    while True:
        line=docu.readline()
        if not line:
            break
        if line=='[Term]\n':
            blocks+=1
            Block=True
            Obsolete=False
            # reset the temporary dictionary
            localDictionary={'is_a':[],'part_of':[],'regulates':[],'positively_regulates':[],'negatively_regulates':[]}
        else :
            if line=='\n':
                Block=False
                # perform a flush of all the informations
                if blocks>1 and not Obsolete:
                    for key in localDictionary.keys():
                        if localDictionary[key]==[]:
                            del localDictionary[key]
                    GO_Terms[localDictionary['id']]=localDictionary
            else:
                if Block:
                    header=line.split(': ')[0].strip()
                    payload=''
                    payload=line.split(': ')[1].strip()
                    if 'CHEBI:' in payload:
                        localDictionary['CHEBI']=payload.split('CHEBI:')[1].split(',')[0].split(']')[0]
                    if header=='id':
                        localDictionary[header]=payload.split(':')[1]
                    if header in ['name', 'namespace', 'def']:
                        localDictionary[header]=payload
                    if header=='is_a':
                        payload=payload.split('!')[0].split(':')[1].strip()
                        localDictionary[header].append(payload)
                    if header=='is_obsolete': 
                        Obsolete=True
                    if header=='relationship':
                        header=str(payload.split()[0].strip())
                        payload=str(payload.split()[1].strip().split(':')[1])
                        localDictionary[header].append(payload)
    
fill_GO_Terms()

