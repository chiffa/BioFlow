"""
The module responsible for parsing of the Uniprot Dataset

Once uniprot is parsed, it is returned as teh dictionary containing the following elements:

Uniprot = { SWISSPROT_ID:{
    'Acnum':[],
    'Names': {'Full': '', 'AltNames': []},
    'GeneRefs': {'Names': [], 'OrderedLocusNames': [], 'ORFNames': []},
    'TaxID': '',
    'Ensembl': [],
    'KEGG': [],
    'EMBL': [],
    'GO': [],
    'Pfam': [],
    'SUPFAM': [],
    'PDB': [],
    'GeneID': [], }}
"""
import BioFlow.configs2 as conf
from BioFlow.Utils.LogManager import logger
import copy

interesting_lines = ['ID', 'AC', 'DE', 'GN', 'OX', 'DR']
interesting_xrefs = ['EMBL', 'GO', 'Pfam', 'Ensembl', 'KEGG', 'PDB', 'GeneID']
names_to_ignore = ['Contains', 'Allergen', 'EC=', 'Flags: ', 'CD_antigen', 'INN=']
starting_dict = {'Acnum': [], 'Names': {'Full': '', 'AltNames': []},
                 'GeneRefs': {'Names': [], 'OrderedLocusNames': [], 'ORFNames': []},
                 'Ensembl': [], 'KEGG': [], 'EMBL': [], 'GO': [], 'Pfam': [], 'SUPFAM': [],
                 'PDB': [], 'GeneID': []}


_ignore = [False, 2]  # a steady constant that regulates a behavior specific to the line skipping


def parse_xref(dico, line):
    """
    Parses an xref line from the Uniprot text file and updates the provided dictionary with the
    results of parsing

    :param dico:
    :param line:
    """
    if 'EMBL; ' in line and 'ChEMBL' not in line:
        contents_list = line.split(';')
        if len(contents_list) > 4:
            package = {'Accession': contents_list[1].strip(),
                       'ID': contents_list[2].strip(),
                       'status': contents_list[3].strip(),
                       'type': contents_list[4].strip().strip('.')}
        else:
            package = {'Accession': contents_list[1].strip(),
                       'ID': contents_list[2].strip(),
                       'status': contents_list[3].strip(),
                       'type': ''}

        dico['EMBL'].append(package)
    if 'GO; GO:' in line:
        dico['GO'].append(line.split(';')[1].split(':')[1].strip())
    if 'Pfam; ' in line:
        dico['Pfam'].append(line.split(';')[1].strip())
    if 'SUPFAM; ' in line:
        dico['SUPFAM'].append(line.split(';')[1].strip())
    if 'Ensembl; ' in line:
        dico['Ensembl'].append(line.split(';')[1].strip())
        dico['Ensembl'].append(line.split(';')[2].strip())
        dico['Ensembl'].append(line.split(';')[3].strip().strip('.'))
    if 'KEGG; ' in line:
        dico['KEGG'].append(line.split(';')[1].strip())
    if 'PDB; ' in line:
        dico['PDB'].append(line.split(';')[1].strip())
    if 'GeneID; ' in line:
        dico['GeneID'].append(line.split(';')[1].strip())


def parse_gene_references(dico, line):
    """
    Parses gene names and references from the UNIPROT text file

    :param dico:
    :param line:
    """
    words = filter(lambda x: x != '', str(line[2:].strip() + ' ').split('; '))
    for word in words:
        if 'ORFNames' in word:
            for subword in word.split('=')[1].strip().split(','):
                dico['GeneRefs']['ORFNames'].append(subword.strip())
        if 'OrderedLocusNames' in word:
            for subword in word.split('=')[1].strip().split(','):
                dico['GeneRefs']['OrderedLocusNames'].append(subword.strip())
        if 'Name=' in word or 'Synonyms=' in word:
            for subword in word.split('=')[1].strip().split(','):
                dico['GeneRefs']['Names'].append(subword.strip())


def parse_name(dico, line):
    """
    Parses a line that contains a name associated to the entry we are trying to load

    :param dico:
    :param line:
    :return:
    """
    if 'RecName: Full=' in line:
        dico['Names']['Full'] = line.split('RecName: Full=')[1].split(';')[0]
        return ''
    if 'AltName: Full=' in line:
        dico['Names']['AltNames'].append(line.split('AltName: Full=')[1].split(';')[0])
        return ''
    if 'Short=' in line:
        dico['Names']['AltNames'].append(line.split('Short=')[1].split(';')[0])
        return ''
    if _ignore[0]:
        if _ignore[1] == 0:
            _ignore[0] = False
            _ignore[1] = 2
            return ''
        else:
            return ''
    if ' Includes:' in line:
        _ignore[0] = True
        return ''
    if any(x in line for x in names_to_ignore):
        return ''         


def process_line(dico, line, keyword):
    """
    A function that processes a line parsed from the UNIPROT database file

    :param dico:
    :param line:
    :param keyword:
    """
    if keyword == 'ID':
        words = filter(lambda a: a != '', line.split(' '))
        dico['ID'] = words[1]
    if keyword == 'AC':
        words = filter(lambda a: a != '', line.split(' '))
        for word in words[1:]:
            dico['Acnum'].append(word.split(';')[0])
    if keyword == 'OX':
        dico['TaxID'] = line.split('NCBI_TaxID=')[1].split(';')[0]
    if keyword == 'DE':
        parse_name(dico, line)
    if keyword == 'GN':
        parse_gene_references(dico, line)
    if keyword == 'DR' and any(x in line for x in interesting_xrefs):
        parse_xref(dico, line)


def end_block(dico, uniprot, tax_id_list):
    """
    Manages the behavior of the end of a parse block

    :param dico:
    :param uniprot:
    :param tax_id_list: list of NCBI taxonomy identifiers wer are interested in
    :return:
    """
    if dico['TaxID'] in tax_id_list:
        _ignore[0] = False
        uniprot[dico['ID']] = dico
    return copy.deepcopy(starting_dict)


def parse_uniprot(source_path=conf.UNIPROT_source, tax_id_to_parse=conf.up_tax_ids):
    """
    Performs the entire uniprot file parsing and importing

    :param source_path: path towards the uniprot test file
    :param tax_id_to_parse: list of NCBI taxonomy identifiers we are interested in
    :return: uniprot parse dictionary
    """
    uniprot = {}
    local_dictionary = copy.deepcopy(starting_dict)
    source_file = open(source_path, "r")
    line_counter = 0
    while True:
        line = source_file.readline()
        line_counter += 1
        if not line:
            break
        keyword = line[0:2]
        if keyword == '//':
            local_dictionary = end_block(local_dictionary, uniprot, tax_id_to_parse)
        if keyword in interesting_lines:
            process_line(local_dictionary, line, keyword)

    logger.info("%s lines scanned during UNIPROT import" % line_counter)
    return uniprot


def get_access_dicts(source_path=conf.UNIPROT_source, tax_id_to_parse=conf.up_tax_ids):
    """
    Returns an access dictionary that would plot genes names, AcNums or EMBL identifiers to the 
    Swissprot IDs

    :param source_path: path towards the uniprot test file
    :param tax_id_to_parse: list of NCBI taxonomy identifiers we are interested in
    :return: dictionary mapping all teh external database identifiers towards uniprot IDs
    """
    uniprot = parse_uniprot(source_path, tax_id_to_parse)
    access_dict = {}
    for key in uniprot.keys():
        for sub_element in uniprot[key]['KEGG']:
            access_dict[sub_element] = key
        for sub_element in uniprot[key]['Ensembl']:
            access_dict[sub_element] = key
        for sub_element in uniprot[key]['EMBL']:
            access_dict[sub_element['Accession']] = key
            access_dict[sub_element['ID']] = key
        for sub_element in uniprot[key]['Acnum']:
            access_dict[sub_element] = key
        for sub_element in uniprot[key]['GeneRefs']['Names']:
            access_dict[sub_element] = key
        for sub_element in uniprot[key]['GeneRefs']['OrderedLocusNames']:
            access_dict[sub_element] = key
        for sub_element in uniprot[key]['GeneRefs']['ORFNames']:
            access_dict[sub_element] = key
    return access_dict


if __name__ == '__main__':
    uniprot = parse_uniprot()
    print len(uniprot)
    print uniprot.keys()
