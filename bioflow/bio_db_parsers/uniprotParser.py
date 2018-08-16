"""
The module responsible for parsing of the Uniprot SWISSPROT .dat file for a subset of
cross-references that are useful in our database.

Once uniprot is parsed, it is returned as the dictionary containing the following elements:

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
import re
import copy
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)

interesting_lines = ['ID', 'AC', 'DE', 'GN', 'OX', 'DR']

interesting_xrefs = ['EMBL', 'GO', 'Pfam', 'Ensembl', 'KEGG', 'PDB', 'GeneID', 'SUPFAM']

names_to_ignore = [
    'Contains',
    'Allergen',
    'EC=',
    'Flags: ',
    'CD_antigen',
    'INN=']

uniprot_load_dict = {
    'Acnum': [],
    'Names': {
        'Full': '',
        'AltNames': []},
    'GeneRefs': {
        'Names': [],
        'AltNames': [],
        'OrderedLocusNames': [],
        'ORFNames': []},
    'Ensembl': [],
    'KEGG': [],
    'EMBL': [],
    'GO': [],
    'Pfam': [],
    'SUPFAM': [],
    'PDB': [],
    'GeneID': [],
    'RefSeq': [],
    'MGI': []}


class UniProtParser(object):
    """Wraps the Uniprot parser """

    def __init__(self, tax_ids_to_parse):
        """

        :param tax_ids_to_parse: list of NCBI taxonomy identifiers we are interested in
        :return:
        """

        self._ignore = [False, 2]
        self.interesting_lines = interesting_lines
        self.interesting_xrefs = interesting_xrefs
        self.names_to_ignore = names_to_ignore
        self._single_up_dict = {}
        self.uniprot = {}
        self.parsed = False
        self.tax_id_list = tax_ids_to_parse

    def parse_xref(self, line):
        """
        Parses an xref line from the Uniprot text file and updates the provided dictionary with the
        results of parsing

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

            self._single_up_dict['EMBL'].append(package)
        if 'GO; GO:' in line:
            self._single_up_dict['GO'].append(line.split(';')[1].split(':')[1].strip())
        if 'Pfam; ' in line:
            self._single_up_dict['Pfam'].append(line.split(';')[1].strip())
        if 'SUPFAM; ' in line:
            self._single_up_dict['SUPFAM'].append(line.split(';')[1].strip())
        if 'Ensembl; ' in line:
            self._single_up_dict['Ensembl'].append(line.split(';')[1].strip())
            self._single_up_dict['Ensembl'].append(line.split(';')[2].strip())
            self._single_up_dict['Ensembl'].append(line.split(';')[3].strip().strip('.'))
        if 'KEGG; ' in line:
            self._single_up_dict['KEGG'].append(line.split(';')[1].strip())
        if 'PDB; ' in line:
            self._single_up_dict['PDB'].append(line.split(';')[1].strip())
        if 'GeneID; ' in line:
            self._single_up_dict['GeneID'].append(line.split(';')[1].strip())
        if 'RefSeq; ' in line:
            self._single_up_dict['RefSeq'].append(line.split(';')[1].strip())
            self._single_up_dict['RefSeq'].append(line.split(';')[2].split(' ')[0].strip())
        if 'MGI;' in line:
            self._single_up_dict['MGI'].append(line.split(';')[2].split(' ')[0].strip())

    def parse_gene_references(self, line):
        """
        Parses gene names and references from the UNIPROT text file

        :param line:
        """
        words = [x for x in str(line[2:].strip() + ' ').split('; ') if x != '']
        for word in words:
            if 'ORFNames' in word:
                for subword in word.split('=')[1].strip().split(','):
                    self._single_up_dict['GeneRefs']['ORFNames'].append(subword.strip())
            if 'OrderedLocusNames' in word:
                for subword in word.split('=')[1].strip().split(','):
                    self._single_up_dict['GeneRefs']['OrderedLocusNames'].append(subword.strip())
            if 'Name=' in word:
                for subword in word.split('=')[1].strip().replace(',', ' ').replace(';', ' ').split():
                    if re.match("^[a-zA-Z0-9_.-]*$", subword):
                        self._single_up_dict['GeneRefs']['Names'].append(subword.strip())
                    else:
                        if '{' not in subword:
                            print "rejected %s: doesn't look like a valid name" % subword
            if 'Synonyms=' in word:
                for subword in word.split('=')[1].strip().replace(',', ' ').replace(';', ' ').split():
                    if re.match("^[a-zA-Z0-9_.-]*$", subword):
                        self._single_up_dict['GeneRefs']['AltNames'].append(subword.strip())
                    else:
                        if '{' not in subword:
                            print "rejected %s: doesn't look like a valid name" % subword

    def parse_name(self, line):
        """
        Parses a line that contains a name associated to the entry we are trying to load

        :param line:
        :return:
        """
        if 'RecName: Full=' in line:
            self._single_up_dict['Names']['Full'] = line.split('RecName: Full=')[1].split(';')[0].split('{')[0]
            return ''
        if 'AltName: Full=' in line:
            self._single_up_dict['Names']['AltNames'].append(
                line.split('AltName: Full=')[1].split(';')[0].split('{')[0])
            return ''
        if 'Short=' in line:
            self._single_up_dict['Names']['AltNames'].append(line.split('Short=')[1].split(';')[0].split('{')[0])
            return ''
        if self._ignore[0]:
            if self._ignore[1] == 0:
                self._ignore[0] = False
                self._ignore[1] = 2
                return ''
            else:
                return ''
        if ' Includes:' in line:
            self._ignore[0] = True
            return ''
        if any(x in line for x in self.names_to_ignore):
            return ''

    def process_line(self, line, keyword):
        """
        A function that processes a line parsed from the UNIPROT database file

        :param line:
        :param keyword:
        """
        if keyword == 'ID':
            words = [a for a in line.split(' ') if a != '']
            self._single_up_dict['ID'] = words[1]
        if keyword == 'AC':
            words = [a for a in line[5:].split(' ') if a != '']
            for word in words:
                self._single_up_dict['Acnum'].append(word.split(';')[0])
        if keyword == 'OX':
            tentative_tax_id = line.split('NCBI_TaxID=')[1].split(';')[0]
            if ' ' in tentative_tax_id:
                tentative_tax_id = tentative_tax_id.split(' ')[0]
            self._single_up_dict['TaxID'] = tentative_tax_id
        if keyword == 'DE':
            self.parse_name(line)
        if keyword == 'GN':
            self.parse_gene_references(line)
        if keyword == 'DR' and any(x in line for x in self.interesting_xrefs):
            self.parse_xref(line)

    def end_block(self):
        """
        Manages the behavior of the end of a parse block

        :return:
        """
        if self._single_up_dict['TaxID'] in self.tax_id_list:
            self._ignore[0] = False
            self.uniprot[self._single_up_dict['ID']] = self._single_up_dict
        return copy.deepcopy(uniprot_load_dict)

    def parse_uniprot(self, source_path):
        """
        Performs the entire uniprot file parsing and importing

        :param source_path: path towards the uniprot test file
        :return: uniprot parse dictionary
        """
        self._single_up_dict = copy.deepcopy(uniprot_load_dict)
        source_file = open(source_path, "r")
        line_counter = 0
        while True:
            line = source_file.readline()
            line_counter += 1
            if not line:
                break
            keyword = line[0:2]
            if keyword == '//':
                self._single_up_dict = self.end_block()
            if keyword in self.interesting_lines:
                self.process_line(line, keyword)

        log.info("%s lines scanned during UNIPROT import", line_counter)
        self.parsed = True
        return self.uniprot

    def get_access_dicts(self):
        """
        Returns an access dictionary that would plot genes names, AcNums or EMBL identifiers to the
        Swissprot IDs

        :return: dictionary mapping all teh external database identifiers towards uniprot IDs
        """
        if not self.parsed:
            log.warning('Attempting to get access points to a non-parsed uniprot object')
        access_dict = {}

        for key in self.uniprot.keys():
            for sub_element in self.uniprot[key]['KEGG']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['Ensembl']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['EMBL']:
                access_dict[sub_element['Accession']] = key
                access_dict[sub_element['ID']] = key
            for sub_element in self.uniprot[key]['Acnum']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['GeneRefs']['Names']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['GeneRefs']['AltNames']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['GeneRefs']['OrderedLocusNames']:
                access_dict[sub_element] = key
            for sub_element in self.uniprot[key]['GeneRefs']['ORFNames']:
                access_dict[sub_element] = key

        return access_dict
