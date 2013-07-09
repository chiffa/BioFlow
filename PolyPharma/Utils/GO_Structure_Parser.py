'''
Created on Jul 2, 2013

@author: andrei
'''

import configs as conf

GO_Terms={}
    # {id:
    # {'id':''
    # 'name':''
    # 'def':''
    # 'namespace':''}}
    #
    # Ignored
    # ['subset', => ignore
    # 'comment', => ignore
    # 'exact_synonym',=> Ignore
    # 'consider', => ignore
    # 'related_synonym', => ignore
    # 'narrow_synonym', => ignore
    # 'broad_synonym', => ignore
    # 'replaced_by', => ignore
    # 'alt_id', => ignore
    # 'xref_analog'] => ignore

GO_Terms_Structure={} 
    # 'is_a'
    # 'relationship',
    #    'part_of'
    #    'regulates'
    #    'positively_regulates'
    #    'negatively_regulates'

'''
Parse dictionary is of the type GO_Term_ID: {}
'''

def fill_GO_Terms():
    docu=open(conf.GeneOntology,"r")
    localDictionary={}
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
            localDictionary={}
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
                        GO_Terms_Structure[localDictionary['id']]=(header,payload)
                    if header=='is_obsolete': 
                        Obsolete=True
                    if header=='relationship':
                        header=str(payload.split()[0].strip())
                        payload=str(payload.split()[1].strip().split(':')[1])
                        GO_Terms_Structure[localDictionary['id']]=(header,payload)
    
fill_GO_Terms()

