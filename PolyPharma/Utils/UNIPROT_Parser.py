'''
Created on Jun 26, 2013

@author: andrei
'''

import configs as conf
import copy 

Interesting_TaxIDs=['9606'] #['36329','9606','1773'] #PLAFA taxonomy id is actually 5833

docu=open(conf.UNIPROT_text,"r")

Interesting_lines=['ID','AC','DE','GN','OX','DR']
Interesing_xrefs=['EMBL','GO','Pfam']

Uniprot={}  # {SWISSPROT_ID:{
            #                 'Acnum':[],
            #                 'RecName':RecordName,
            #                 'ORFNames': [],
            #                 'OrderedLocusNames':[],
            #                 'TaxID': NCBI_TaxID,
            #                 'EMBL':[],
            #                 'GO':[],
            #                 'Pfam':[],
            #                 }}

defDict={'Flags':[],'Acnum':[],'Names':{'Full':'','AltNames':[],'Includes':[]},'GeneRefs':{'Names':[],'OrderedLocusNames':[],'ORFNames':[]},'EMBL':[],'GO':[],'Pfam':[], 'SUPFAM':[]}

Buffer_vars=[False,[]]

def parse_Xref(Dico,Line):
    if 'EMBL; ' in Line and 'ChEMBL' not in Line:
        splt=Line.split(';')
        package=(splt[1].strip(),splt[2].strip(),splt[3].strip())
        Dico['EMBL'].append(package)
    if 'GO; GO:' in Line:
        Dico['GO'].append(Line.split(';')[1].split(':')[1].strip())
    if 'Pfam; ' in Line:
        Dico['Pfam'].append(Line.split(';')[1].strip())
    if 'SUPFAM; ' in Line:
        Dico['SUPFAM'].append(Line.split(';')[1].strip())

def parse_GeneRefs(Dico,Line):
    words=filter(lambda a:a!='', str(Line.strip()+' ').split('; '))
    for word in words[1:]:
        if 'ORFNames' in word:
            for subword in word.split('=')[1].strip().split(','):
                Dico['GeneRefs']['ORFNames'].append(subword.strip())
        if 'OrderedLocusNames' in word:
            for subword in word.split('=')[1].strip().split(','):
                Dico['GeneRefs']['OrderedLocusNames'].append(subword.strip())
        if 'Name=' in word or 'Synonyms=' in word:
            for subword in word.split('=')[1].strip().split(','):
                Dico['GeneRefs']['Names'].append(subword.strip())


def parse_Name(Dico,Line):
    # TODO: improve the includes and contains protections of variables => We don't need them so far, but it is still better
    # TODO: add FlagIgnore
    if 'RecName: Full=' in Line:
        if Buffer_vars[0]==False:
            Dico['Names']['Full']=Line.split('RecName: Full=')[1].split(';')[0]
        else:
            Buffer_vars[1].append(Line.split('RecName: Full=')[1].split(';')[0])
    else:
        if 'AltName: Full=' in Line:
            Dico['Names']['AltNames'].append(Line.split('AltName: Full=')[1].split(';')[0])
        else:
            if 'Short=' in Line:
                Dico['Names']['AltNames'].append(Line.split('Short=')[1].split(';')[0])
            else:
                if 'EC=' in Line:
                    Buffer_vars[1].append(Line.split('EC=')[1].split(';')[0])
                else:
                    if ' Includes:' in Line:
                        Buffer_vars[0]=True
                        if Buffer_vars[1]!=[]:
                            Dico['Names']['Includes'].append(Buffer_vars[1])
                            Buffer_vars[1]=[]
#                     else:
#                         print 'issue1 here \n',Line,'\n','<========================>'
                
                
def process_line(Dico, Line, keyword):
    if keyword=='ID':
        words=filter(lambda a:a!='', Line.split(' '))
        Dico['ID']=words[1]
    if keyword=='AC':
        words=filter(lambda a:a!='', Line.split(' '))
        for word in words[1:]:
            Dico['Acnum'].append(word)
    if keyword=='OX':
        Dico['TaxID']=Line.split('NCBI_TaxID=')[1].split(';')[0]
    if keyword=='DE':
        parse_Name(Dico,Line)
    if keyword=='GN':
        parse_GeneRefs(Dico,Line)
    if keyword=='DR' and any(x in Line for x in Interesing_xrefs):
        parse_Xref(Dico,Line)


def end_Block(Dico):
    if Dico['TaxID'] in Interesting_TaxIDs:
        Buffer_vars[0]=False
        Uniprot[Dico['ID']]=Dico
    return copy.deepcopy(defDict)

def Parse_Uniprot():
    LocalDictionary=copy.deepcopy(defDict)
    while True:
        line = docu.readline()
        if not line:
            break
        keyword=line[0:2]
        if keyword=='//':
            LocalDictionary=end_Block(LocalDictionary)
        if  keyword in Interesting_lines:
            process_line(LocalDictionary, line, keyword)

Parse_Uniprot()
print len(Uniprot)
