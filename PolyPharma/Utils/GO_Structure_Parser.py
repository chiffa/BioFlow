'''
Created on Jul 2, 2013

@author: andrei
'''

import PolyPharma.configs as conf


def fill_GO_Terms():
    """
    Parse dictionary is of the type GO_Term_ID: {}
    """
    GO_Terms = {}
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

    GO_Terms_Structure = [] # [(Node1,relation,Node2)]
        #  Where relation in:
        # 'is_a'
        # 'relationship',
        #    'part_of'
        #    'regulates'
        #    'positively_regulates'
        #    'negatively_regulates'

    docu = open(conf.GeneOntology, "r")
    localDictionary = {}
    LocalRelations = []
    blocks = 0
    Block = False
    Obsolete = False
    while True:
        line = docu.readline()
        if not line:
            break
        if line == '[Term]\n':
            blocks += 1
            Block = True
            Obsolete = False
            # reset the temporary dictionary
            localDictionary = {}
            LocalRelations = []
        else :
            if line == '\n':
                Block = False
                # perform a flush of all the informations
                if blocks > 1 and not Obsolete:
                    for key in localDictionary.keys():
                        if localDictionary[key] == []:
                            del localDictionary[key]
                    GO_Terms[localDictionary['id']] = localDictionary
                    GO_Terms_Structure = GO_Terms_Structure + LocalRelations
            else:
                if Block:
                    header = line.split(': ')[0].strip()
                    payload = ''
                    payload = line.split(': ')[1].strip()
                    if 'CHEBI:' in payload:
                        localDictionary['CHEBI'] = payload.split('CHEBI:')[1].split(',')[0].split(']')[0]
                    if header == 'id':
                        localDictionary[header] = payload.split(':')[1]
                    if header in ['name', 'namespace', 'def']:
                        localDictionary[header] = payload
                    if header == 'is_a':
                        payload = payload.split('!')[0].split(':')[1].strip()
                        LocalRelations.append((localDictionary['id'], header, payload))
                    if header == 'is_obsolete':
                        Obsolete = True
                    if header == 'relationship':
                        header = str(payload.split()[0].strip())
                        payload = str(payload.split()[1].strip().split(':')[1])
                        LocalRelations.append((localDictionary['id'], header, payload))
    return GO_Terms, GO_Terms_Structure

