"""
Module containing the Reactome Biopax lvl3 .owl file parser
"""
import xml.etree.ElementTree as ET
from collections import defaultdict
from bioflow.utils.log_behavior import get_logger
from bioflow.configs import main_configs
import os


log = get_logger(__name__)


def zip_dicts(dict1, dict2):
    """
    defines a dictionary update rules that are required for a proper dictionary update

    :param dict1: dictionary #1
    :param dict2: dictionary #2
    """
    for key in list(dict2.keys()):
        if key not in list(dict1.keys()):
            dict1[key] = dict2[key]  # never used in production
        else:
            assert isinstance(dict2[key], (list, tuple))
            dict1[key] = dict1[key] + dict2[key]
    return dict1


class ReactomeParser(object):
    """
    Wrapper class for the Reactome parser routines

    Documentation for the contents that are filled by a parse:

        - self.BioSources = {ID:name}
        - self.CellularLocations = {Id:name}
        - self.SeqModVoc = {Id:name}
        - self.SeqSite = {Id:Position}
        - self.Pathways = {ID:{'displayName':'','pathwayComponent':'','name':[], 'PathwayStep':[]}}
        - self.PathwaySteps = {ID:{'stepProcess':[],'nextStep':[]}}

        - self.DnaRefs = {ID:{'name':[],'ENSEMBL':'','organism':''}}
        - self.RnaRefs = {ID:{'name':[],'ENSEMBL':'', 'miRBase':'', EMBL:'', s'organism':''}}
        - self.SmallMoleculeRefs = {ID:{'name':[],'ChEBI':'','organism':''}}
        - self.ProteinRefs = {ID:{'name':[],'UniProt':'','organism':''}}
        - self.ModificationFeatures = {ID:{'location':'','modification':''}}

        - self.TemplateReactions = {ID:{'product':'','displayName':'',
                                    'references':{'names':[], ...}}}
        - self.Degradations = {ID:{'product':'','displayName':'',
                                'references':{'eCNumber':[],...}}}
        - self.BiochemicalReactions = {ID:{'left':[],'right':[],'displayName':'',
                                        'references':{'eCNumber':[],...}}}
        - self.Catalysises = {ID:{Controller:'', Controlled:'', controlType:''}}
        - self.Modulations ={ID:{modulator, modulated}
        ( This is essentially a compressed regulation of activity of the catalysts)

    All the elements in the physical entities and physical entities collection follow the
    same pattern:
    {Id:{'cellularLocation':'', displayName:'', collectionMembers':[], 'references':{
    'names':[],...}}}

    Those entities include:

        self.Dnas
        self.Dna_Collections
        self.Rnas
        self.Rna_Collections
        self.SmallMolecules
        self.SmallMolecule_Collections
        self.Proteins
        self.Protein_Collections
        self.PhysicalEntities
        self.PhysicalEntity_Collections
        self.Complexes
        self.Complex_Collections
    """

    # TODO: modify external elements to remove legacy methods support
    # TODO: add unification x-ref parsing to match things like cellular locations to GO terms
    # TODO: move the definition of the file to parse from the __inti__ method to the execution tree
    # TODO: add the following post-processing filters: , ['BiochemicalReaction'] on Catalysis
    #  when parsing 'Control' object set
    # REFACTOR: [MAINTENABILITY]: convert the system of references to a single tree that has all
    #  the different subcategories already declared

    def __init__(self, path_to_biopax_file=main_configs.reactome_biopax_path):

        self.tree = ET.parse(path_to_biopax_file)
        self.root = self.tree.getroot()
        log.info('Reactome parser parsed the xml tree')

        self.BioSources = {}
        self.CellularLocations = {}
        self.SeqModVoc = {}
        self.SeqSite = {}
        self.Pathways = {}
        self.PathwaySteps = {}

        self.DnaRefs = {}
        self.RnaRefs = {}
        self.SmallMoleculeRefs = {}
        self.ProteinRefs = {}
        self.ModificationFeatures = {}

        self.Dnas = {}
        self.Dna_Collections = {}
        self.Rnas = {}
        self.Rna_Collections = {}
        self.SmallMolecules = {}
        self.SmallMolecule_Collections = {}
        self.Proteins = {}
        self.Protein_Collections = {}
        self.PhysicalEntities = {}
        self.PhysicalEntity_Collections = {}
        self.Complexes = {}
        self.Complex_Collections = {}

        self.TemplateReactions = {}
        self.Degradations = {}
        self.BiochemicalReactions = {}
        self.Catalysises = {}
        self.Modulations = {}

        self.parsed = False

    def _find_in_root(self, term_name):
        """
        Abstracts the simplest xml foot finding of a path. Currently a test function to be
        integrated in more high-level functions

        :param term_name:
        :return: an xml iterator
        """
        search_pattern = '{http://www.biopax.org/release/biopax-level3.owl#}%s' % term_name
        xml_iterator = self.root.findall(search_pattern)
        return xml_iterator

    def _single_tag_parse(self, primary_term, target_dict, tag_to_parse):
        """
        Parses a simple base objects in the Reactome lvl3 Biopax .owl file

        Used for the most basic parsing routines within the Reactome

        :param primary_term: term we want to anchor to for parsing
        :param target_dict: dict we want to insert the parse results into
        :param tag_to_parse: what tag we want to parse first
        """
        for object_of_interest in self._find_in_root(primary_term):
            key_ = list(object_of_interest.attrib.values())[0]
            for object_property in object_of_interest:
                if tag_to_parse in object_property.tag:
                    target_dict[key_] = object_property.text

    def _parse_xref(self, primary_term, target_dict, mapped_terms):
        """
        Parses xref type objects in the Reactome lvl3 Biopax .owl file

        :param primary_term: term we want to anchor to for parsing
        :param target_dict: dict we want to insert the parse results into
        :param mapped_terms: specific x-refs whose content we would want to extract
        """
        for object_of_interest in self._find_in_root(primary_term):
            key_ = list(object_of_interest.attrib.values())[0]
            target_dict[key_] = {'name': []}
            for object_property in object_of_interest:
                if '}name' in object_property.tag and object_property.text:
                    other_term_inserted = False
                    for term in mapped_terms:
                        if term in object_property.text:
                            for word in object_property.text.split():
                                if term in word:
                                    tag = word.split(':')[0].strip('[]')
                                    tagged = word.split(':')[1].strip('[]')
                                    target_dict[key_][tag] = tagged
                                    other_term_inserted = True
                    if not other_term_inserted:
                        target_dict[key_]['name'].append(object_property.text)
                if '}organism' in object_property.tag \
                        and list(object_property.attrib.values())[0][1:] != 'BioSource1':
                    # second condition removes references to the main organism in the xml file
                    target_dict[key_]['organism'] = self.BioSources[
                        list(object_property.attrib.values())[0][1:]]
            target_dict[key_]['name'] = list(set(target_dict[key_]['name']))

    def _parse_modification_features(self):
        """
        Parses modification features.

        Because of unicity of parsing syntax they cannot be folded elsewhere
        """
        for single_ModificationFeature in self._find_in_root('ModificationFeature'):
            key = list(single_ModificationFeature.attrib.values())[0]
            self.ModificationFeatures[key] = {'ID': key}  # pre-processing
            for modification_property in single_ModificationFeature:
                if '}featureLocation' in modification_property.tag:  # trigger #1
                    # dict_position = element to retrieve
                    self.ModificationFeatures[key]['location'] = self.SeqSite[
                        list(modification_property.attrib.values())[0][1:]]  # trigger #2
                if '}modificationType' in modification_property.tag:
                    # dict_position = element to retrieve
                    self.ModificationFeatures[key]['modification'] = self.SeqModVoc[
                        list(modification_property.attrib.values())[0][1:]]

    def _reactome_obj_parse_tag(self, local_dict, local_property, reactome_xrefs, is_collection):
        """
        Parses a meta-physical entity tag.

        :param local_dict: local dictionary where the parser needs to be returned
        :param local_property: property we are currently processing
        :param reactome_xrefs: dictionary of x-references, for reference tag parsing
        :param is_collection: Fall-through tag indicating if we are parsing a collection
        """
        if '}cellularLocation' in local_property.tag:
            local_dict['cellularLocation'] = list(local_property.attrib.values())[0][1:]
        if '}displayName' in local_property.tag:
            local_dict['displayName'] = local_property.text
        if '}name' in local_property.tag:
            local_dict['references']['name'].append(local_property.text)
        if '}memberPhysicalEntity' in local_property.tag:
            is_collection = True
            local_dict['collectionMembers'].append(
                list(local_property.attrib.values())[0][1:])
        if '}entityReference' in local_property.tag:
            # print local_property.attrib.values()[0][1:]
            # TODO: temporary fix < HAHAHA
            if not "EntityReference" in list(local_property.attrib.values())[0][1:]:
                # pass
                # print local_property.attrib.values()[0][1:]
                # for ref in reactome_xrefs.keys():
                #     print ref
                # meta_ref_types = set(''.join([i for i in ref if not i.isdigit()]) for ref in reactome_xrefs.keys())
                # print meta_ref_types
                local_dict['references'] = \
                        zip_dicts(reactome_xrefs[list(local_property.attrib.values())[0][1:]],
                                  local_dict['references'])
        if '}feature' in local_property.tag \
                and 'ModificationFeature' in list(local_property.attrib.values())[0]:
            local_dict['modification'].append(self.ModificationFeatures[
                                                  list(local_property.attrib.values())[0][1:]])
        if '}component' in local_property.tag \
                and 'Stoichiometry' not in local_property.tag:
                    local_dict['parts'].append(
                        list(local_property.attrib.values())[0][1:])
        return is_collection

    def _reactome_obj_parse(self, primary_term, reactome_obj_dict,
                            reactome_obj_collection_dict, reactome_obj_xrefs,
                            sup_mods=False, sup_parts=False):
        """
        A common parse round for all the meta-physical entity

        :param primary_term: term name within Reactome file
        :param reactome_obj_dict: dictionary of meta-object parses where we are inserting it
        :param reactome_obj_collection_dict: dictionary of collections of meta-objects
        :param reactome_obj_xrefs: dictionary containing the x-refs for the meta-object
        :param sup_mods: if true, modulation will be suppressed
        :param sup_parts: if true, parts tag will be suppressed
        """
        for reactome_object in self._find_in_root(primary_term):

            key_ = list(reactome_object.attrib.values())[0]
            base_dict = {'collectionMembers': [],
                         'modification': [],
                         'parts': [],
                         'references': {
                             'name': []}
                         }
            is_collection = False

            for reactome_object_property in reactome_object:
                is_collection = self._reactome_obj_parse_tag(base_dict, reactome_object_property,
                                                             reactome_obj_xrefs, is_collection)

            base_dict['references']['name'] = list(set(base_dict['references']['name']))
            base_dict['modification'] = [dict(t)
                                         for t in set([tuple(d.items())
                                                       for d in base_dict['modification']])]

            if is_collection:
                del base_dict['modification']
                del base_dict['parts']
                reactome_obj_collection_dict[key_] = base_dict

            else:
                del base_dict['collectionMembers']
                if not base_dict['modification'] or sup_mods:
                    del base_dict['modification']
                if not base_dict['parts'] or sup_parts:
                    del base_dict['parts']
                reactome_obj_dict[key_] = base_dict

    def _parse_reaction(self, primary_term, target_dict, tags_to_parse=()):
        """
        A meta-parser for the reactions type of data

        :param primary_term: term within Reactome file we would like to match
        :param target_dict: dictionary where we will be dumping the parse results
        :param tags_to_parse: tags for a Reactome object we would like to parse
        :param flatten: if we want to flatten left/right tags parse results (legacy support)
        :param remap: if we want to rename some terms we've parsed (legacy support)
        """
        for reaction_object in self._find_in_root(primary_term):
            key_ = list(reaction_object.attrib.values())[0]
            base_dict = {'right': [],
                         'left': [],
                         'references': {
                             'name': []}}
            for reaction_property in reaction_object:
                if '}product' in reaction_property.tag:
                    base_dict['right'].append(list(reaction_property.attrib.values())[0][1:])
                if '}displayName' in reaction_property.tag:
                    base_dict['displayName'] = reaction_property.text
                if '}name' in reaction_property.tag:
                    base_dict['references'][
                        'name'] = reaction_property.text
                if '}eCNumber' in reaction_property.tag:
                    base_dict['references']['eCNumber'] = reaction_property.text
                if '}left' in reaction_property.tag:
                    base_dict['left'].append(list(reaction_property.attrib.values())[0][1:])
                if '}right' in reaction_property.tag:
                    base_dict['right'].append(list(reaction_property.attrib.values())[0][1:])

            for key in list(base_dict.keys()):
                if key not in tags_to_parse and key not in ['references', 'displayName']:
                    del base_dict[key]

            target_dict[key_] = base_dict

    def _parse_a_catalysis(self, primary_term, target_dict):
        """
        General function of catalysis parsing

        :param primary_term: term in the Reactome we would like to parse
        :param target_dict: dictionary into which we would like to insert parse results
        (legacy reasons)
        """
        for catalysis_object in self._find_in_root(primary_term):
            key_ = list(catalysis_object.attrib.values())[0]
            base_dict = {}
            for catalysis_property in catalysis_object:
                if '}controlled' in catalysis_property.tag:
                    base_dict['controlled'] = list(catalysis_property.attrib.values())[0][1:]
                if '}controller' in catalysis_property.tag:
                    base_dict['controller'] = list(catalysis_property.attrib.values())[0][1:]
                if '}controlType' in catalysis_property.tag:
                    base_dict['ControlType'] = catalysis_property.text

            if 'Pathway' not in base_dict['controlled']:
                target_dict[key_] = base_dict

    def _parse_modulations(self):
        """
        Parses modulations
        """
        for single_Modulation in self._find_in_root('Modulation'):
            key = list(single_Modulation.attrib.values())[0]
            self.Modulations[key] = {}
            for modulation_property in single_Modulation:
                if '}displayName' in modulation_property.tag:
                    self.Modulations[key]['displayName'] = modulation_property.text
                if '}controller' in modulation_property.tag:
                    self.Modulations[key]['controller'] = list(modulation_property.attrib.values())[0][
                        1:]
                if '}controlled' in modulation_property.tag:
                    self.Modulations[key]['controlled'] = \
                        self.Catalysises[list(modulation_property.attrib.values())[0][1:]]['controller']
                if '}controlType' in modulation_property.tag:
                    self.Modulations[key]['controlType'] = modulation_property.text

    def _parse_pathways(self):
        """
        Parses Pathways
        """
        for single_Pathway in self._find_in_root('Pathway'):
            key = list(single_Pathway.attrib.values())[0]
            local_dict = {
                'components': [],
                'references': {
                    'name': []},
                'PathwayStep': []}
            for pathway_property in single_Pathway:
                if '}displayName' in pathway_property.tag:
                    local_dict['displayName'] = pathway_property.text
                if '}name' in pathway_property.tag:
                    local_dict['references']['name'].append(pathway_property.text)
                if '}pathwayComponent' in pathway_property.tag:
                        # and 'Pathway' in list(pathway_property.attrib.values())[0]:
                        # only will parse things that has "Pathway in them => disabled
                    local_dict['components'].append(
                        list(pathway_property.attrib.values())[0][1:])
                if '}pathwayOrder' in pathway_property.tag:
                    local_dict['PathwayStep'].append(
                        list(pathway_property.attrib.values())[0][1:])
            self.Pathways[key] = local_dict

    def _parse_pathway_steps(self):
        """
        Parses Pathway steps
        """
        exclude = [
            'Modulation',  # does not exist anymore
            'Control',
            'TemplateReactionRegulation',  # does not exist anymore
            'Catalysis']

        for single_Pathway_step in self._find_in_root('PathwayStep'):
            key = list(single_Pathway_step.attrib.values())[0]
            local_dict = {'components': [], 'nextStep': []}
            for pathway_property in single_Pathway_step:
                if '}stepProcess' in pathway_property.tag and not any(
                        x in list(pathway_property.attrib.values())[0] for x in exclude):
                    local_dict['components'].append(
                        list(pathway_property.attrib.values())[0][1:])
                if '}nextStep' in pathway_property.tag:
                    local_dict['nextStep'].append(
                        list(pathway_property.attrib.values())[0][1:])
            self.PathwaySteps[key] = local_dict

    def parse_all(self):
        """
        The only method that should be called publicly to ensure everything was parsed and returned
        properly
        """
        self._single_tag_parse('BioSource', self.BioSources, '}name')
        self._single_tag_parse('CellularLocationVocabulary', self.CellularLocations, '}term')
        self._single_tag_parse('SequenceModificationVocabulary', self.SeqModVoc, '}term')
        self._single_tag_parse('SequenceSite', self.SeqSite, '}sequencePosition')

        self._parse_xref('DnaReference', self.DnaRefs, ['ENSEMBL:'])
        self._parse_xref('RnaReference', self.RnaRefs, ['ENSEMBL:', 'miRBase:', 'EMBL:'])
        self._parse_xref('SmallMoleculeReference', self.SmallMoleculeRefs, ['ChEBI:'])
        self._parse_xref('ProteinReference', self.ProteinRefs, ['UniProt:'])

        self._parse_modification_features()

        self._reactome_obj_parse('Dna', self.Dnas, self.Dna_Collections, self.DnaRefs,
                                 sup_mods=True, sup_parts=True)
        self._reactome_obj_parse('Rna', self.Rnas, self.Rna_Collections, self.RnaRefs,
                                 sup_mods=True, sup_parts=True)
        self._reactome_obj_parse('SmallMolecule', self.SmallMolecules, self.SmallMolecule_Collections,
                                 self.SmallMoleculeRefs, sup_mods=True, sup_parts=True)
        self._reactome_obj_parse('Protein', self.Proteins, self.Protein_Collections,
                                 self.ProteinRefs, sup_mods=False, sup_parts=True)
        self._reactome_obj_parse('PhysicalEntity', self.PhysicalEntities, self.PhysicalEntity_Collections,
                                 defaultdict(None), sup_mods=True, sup_parts=True)
        self._reactome_obj_parse('Complex', self.Complexes, self.Complex_Collections,
                                 defaultdict(None), sup_mods=True, sup_parts=False)

        self._parse_reaction('TemplateReaction', self.TemplateReactions, ['product'])
        self._parse_reaction('Degradation', self.Degradations, ['left'])
        self._parse_reaction('BiochemicalReaction', self.BiochemicalReactions, ['left', 'right'])

        self._parse_a_catalysis('Catalysis', self.Catalysises)
        # self._parse_a_catalysis('TemplateReactionRegulation', self.Catalysises)  # Does not exist
        # # anymore
        self._parse_a_catalysis('Control', self.Catalysises)

        # self._parse_modulations()  # ceased to exist
        self._parse_pathways()
        self._parse_pathway_steps()

        self.parsed = True

        log.info('Reactome parser finished parsing xml tree to dict collection')


if __name__ == "__main__":
    source_file = 'unittests/UT_examples/reactome_extract.owl'
    # log.debug(os.path.abspath(source_file))
    # print(os.path.abspath(source_file))
    # source_file = main_configs.reactome_biopax_path
    RP = ReactomeParser(source_file)
    RP.parse_all()
