'''
Created on Jun 26, 2013

@author: andrei
'''

import PolyPharma.configs as conf
import copy

Interesting_TaxIDs = [TaxID.strip() for TaxID in conf.Sources['UNIPROT']['tax_ids'].split(',') if TaxID not in ('',' ')]
Interesting_lines = ['ID', 'AC', 'DE', 'GN', 'OX', 'DR']
Interesing_xrefs = ['EMBL', 'GO', 'Pfam', 'Ensembl', 'KEGG', 'PDB', 'GeneID']
NameIgnore = ['Contains', 'Allergen', 'EC=', 'Flags: ', 'CD_antigen', 'INN=']
defDict = {'Acnum':[], 'Names':{'Full':'', 'AltNames':[]}, 'GeneRefs':{'Names':[], 'OrderedLocusNames':[], 'ORFNames':[]},
           'Ensembl':[], 'KEGG':[], 'EMBL':[], 'GO':[], 'Pfam':[], 'SUPFAM':[], 'PDB':[], 'GeneID':[]}


lst1 = [
        'YOR031W',
        'YOR001W',
        'YOL107W',
        'YOL124C',
        'YOL040C',
        'YOR184W',
        'YOR374W',
        'YOR125C',
        'YOL087C',
        'YOR334W',
    ]


# TODO: refactor to avoid any call to the expensive function unless a specific function has been build

# Uniprot = {}  # {SWISSPROT_ID:{
            #                 'Acnum':[],
            #                 'RecName':RecordName,
            #                 'ORFNames': [],
            #                 'OrderedLocusNames':[],
            #                 'TaxID': NCBI_TaxID,
            #                 'EMBL':[],
            #                 'GO':[],
            #                 'Pfam':[],
            #                 }}

Ignore = [False, 2]


def parse_Xref(Dico,Line):
    """

    :param Dico:
    :param Line:
    """
    if 'EMBL; ' in Line and 'ChEMBL' not in Line:
        splt = Line.split(';')
        if len(splt)>4:
            package = { 'Accession': splt[1].strip(), 'ID':splt[2].strip(), 'status':splt[3].strip(), 'type':splt[4].strip().strip('.')}
        else:
            package = { 'Accession': splt[1].strip(), 'ID':splt[2].strip(), 'status':splt[3].strip(), 'type':''}

        Dico['EMBL'].append(package)
    if 'GO; GO:' in Line:
        Dico['GO'].append(Line.split(';')[1].split(':')[1].strip())
    if 'Pfam; ' in Line:
        Dico['Pfam'].append(Line.split(';')[1].strip())
    if 'SUPFAM; ' in Line:
        Dico['SUPFAM'].append(Line.split(';')[1].strip())
    if 'Ensembl; ' in Line:
        Dico['Ensembl'].append(Line.split(';')[1].strip())
        Dico['Ensembl'].append(Line.split(';')[2].strip())
        Dico['Ensembl'].append(Line.split(';')[3].strip().strip('.'))
    if 'KEGG; ' in Line:
        Dico['KEGG'].append(Line.split(';')[1].strip())
    if 'PDB; ' in Line:
        Dico['PDB'].append(Line.split(';')[1].strip())
    if 'GeneID; ' in Line:
        Dico['GeneID'].append(Line.split(';')[1].strip())


def parse_GeneRefs(Dico,Line):
    """

    :param Dico:
    :param Line:
    """
    words = filter(lambda a:a != '', str(Line[2:].strip() + ' ').split('; '))
    for word in words:
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
    """

    :param Dico:
    :param Line:
    :return:
    """
    if 'RecName: Full=' in Line:
        Dico['Names']['Full'] = Line.split('RecName: Full=')[1].split(';')[0]
        return ''
    if 'AltName: Full=' in Line:
        Dico['Names']['AltNames'].append(Line.split('AltName: Full=')[1].split(';')[0])
        return ''
    if 'Short=' in Line:
        Dico['Names']['AltNames'].append(Line.split('Short=')[1].split(';')[0])
        return ''
    if Ignore[0]:
        if Ignore[1] == 0:
            Ignore[0] = False
            Ignore[1] = 2
            return ''
        else:
            return ''
    if ' Includes:' in Line:
        Ignore[0] = True
        return ''
    if any(x in Line for x in NameIgnore):
        return ''         


def process_line(Dico, Line, keyword):
    """

    :param Dico:
    :param Line:
    :param keyword:
    """
    if keyword == 'ID':
        words = filter(lambda a:a != '', Line.split(' '))
        Dico['ID'] = words[1]
    if keyword == 'AC':
        words = filter(lambda a:a != '', Line.split(' '))
        for word in words[1:]:
            Dico['Acnum'].append(word.split(';')[0])
    if keyword == 'OX':
        Dico['TaxID']=Line.split('NCBI_TaxID=')[1].split(';')[0]
    if keyword == 'DE':
        parse_Name(Dico,Line)
    if keyword == 'GN':
        parse_GeneRefs(Dico,Line)
    if keyword == 'DR' and any(x in Line for x in Interesing_xrefs):
        parse_Xref(Dico,Line)


def end_Block(Dico, Uniprot):
    """

    :param Dico:
    :return:
    """
    if Dico['TaxID'] in Interesting_TaxIDs:
        Ignore[0] = False
        Uniprot[Dico['ID']] = Dico
    return copy.deepcopy(defDict)


def Parse_Uniprot():
    """


    """
    Uniprot = {}
    LocalDictionary = copy.deepcopy(defDict)
    source_file = open(conf.UNIPROT_source, "r")
    while True:
        line = source_file.readline()
        if not line:
            break
        keyword = line[0:2]
        if keyword == '//':
            LocalDictionary = end_Block(LocalDictionary, Uniprot)
        if  keyword in Interesting_lines:
            process_line(LocalDictionary, line, keyword)

    return Uniprot


def get_Names_dict():
    """

    :param Uniprot:
    :return:
    """
    Uniprot = Parse_Uniprot()
    namesDict = {}
    for elt in Uniprot.keys():
        NameList = [Uniprot[elt]['Names']['Full'].lower().strip()]
        for subelt in Uniprot[elt]['Names']['AltNames']:
            NameList.append(subelt.lower().strip())
        for subelt in Uniprot[elt]['GeneRefs']['Names']:
            NameList.append(subelt.lower().strip())
        for subelt in Uniprot[elt]['GeneRefs']['OrderedLocusNames']:
            NameList.append(subelt.lower().strip())
        for subelt in Uniprot[elt]['GeneRefs']['ORFNames']:
            NameList.append(subelt.lower().strip())
        
        for subelt in NameList:
            namesDict[subelt] = elt
    
    return namesDict


def get_access_dicts():
    '''
    Returns an access dictionary that would plot genes names, AcNums or EMBL identifiers to the 
    Swissprot IDs

    :param Uniprot:
    :return:
    '''
    Uniprot = Parse_Uniprot()
    access_dict = {}
    for key in Uniprot.keys():
        for subelt in Uniprot[key]['KEGG']:
            access_dict[subelt] = key
        for subelt in Uniprot[key]['Ensembl']:
            access_dict[subelt] = key
        for subelt in Uniprot[key]['EMBL']:
            access_dict[subelt['Accession']] = key
            access_dict[subelt['ID']] = key
        for subelt in Uniprot[key]['Acnum']:
            access_dict[subelt] = key
        for subelt in Uniprot[key]['GeneRefs']['Names']:
            access_dict[subelt] = key
        for subelt in Uniprot[key]['GeneRefs']['OrderedLocusNames']:
            access_dict[subelt] = key
        for subelt in Uniprot[key]['GeneRefs']['ORFNames']:
            access_dict[subelt] = key
    return access_dict

if __name__ == '__main__':
    Uniprot = Parse_Uniprot()
    print len(Uniprot)
    # names_Dict = get_Names_dict()
    # print len(names_Dict)
    # access_dict = get_access_dicts()
    # print len(access_dict)
    # accumulator = []
    # for id, aclist in Uniprot.iteritems():
    #     for item in aclist['Acnum']:
    #         accumulator.append(item)
    # print accumulator