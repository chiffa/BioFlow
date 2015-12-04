"""
Contains the functions responsible for the parsing of the GO terms
"""
from BioFlow.utils.log_behavior import logger


class GOTermsParser(object):
    """ Wrapper object for a parser of GO terms."""

    def __init__(self):
        self.go_terms = {}
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

        self.go_terms_structure = []
        # [(Node1, relation, Node2)]
        #  Where relation in:
        # 'is_a'
        # 'relationship',
        #    'part_of'
        #    'regulates'
        #    'positively_regulates'
        #    'negatively_regulates'

        self.local_dictionary = {}
        self.local_relations = []
        self.blocks = 0
        self.block = False
        self.obsolete = False

    def start_block(self):
        """ resets temporary stores to fill so that a new term can be loaded """
        self.blocks += 1
        self.block = True
        self.obsolete = False
        self.local_dictionary = {}
        self.local_relations = []

    def parse_line_in_block(self, header, payload):
        """
        Parses a line within GO term parameters block

        :param header: GO term parameter name
        :param payload: GO term parameter value
        """
        if 'CHEBI:' in payload:
            self.local_dictionary['CHEBI'] = payload.split(
                'CHEBI:')[1].split(',')[0].split(']')[0]
        if header == 'id':
            self.local_dictionary[header] = payload.split(':')[1]
        if header in ['name', 'namespace', 'def']:
            self.local_dictionary[header] = payload
        if header == 'is_a':
            payload = payload.split('!')[0].split(':')[1].strip()
            self.local_relations.append(
                (self.local_dictionary['id'], header, payload))
        if header == 'is_obsolete':
            self.obsolete = True
        if header == 'relationship':
            header = str(payload.split()[0].strip())
            payload = str(payload.split()[1].strip().split(':')[1])
            self.local_relations.append(
                (self.local_dictionary['id'], header, payload))

    def flush_block(self):
        """ flushes all temporary term stores to the main data stores """
        self.block = False
        if not self.obsolete and self.local_dictionary:
                self.go_terms[self.local_dictionary['id']] = self.local_dictionary
                self.go_terms_structure = self.go_terms_structure + self.local_relations

    def parse_go_terms(self, source_file_path):
        """
        Takes the path to the gene ontology .obo file and returns result of parse dict and list

        :param source_file_path: gene ontology .obo file
        :return: dict containing term parse, list containing inter-term relationship (turtle)
        triplets
        """
        with open(source_file_path, "r") as go_terms_source:
            for line in go_terms_source:
                if line == '[Term]\n':
                    self.start_block()
                elif self.block and line == '\n':
                    self.flush_block()
                elif self.block:
                    try:
                        header = line.split(': ')[0].strip()
                        payload = line.split(': ')[1].strip()
                        self.parse_line_in_block(header, payload)
                    except IndexError:
                        logger.error("Line '%s' violates obo conventions." % line +
                                     " Please check file integrity")

        return self.go_terms, self.go_terms_structure
