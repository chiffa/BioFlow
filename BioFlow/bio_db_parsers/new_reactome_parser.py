"""
Module containing the reactome parser routines.
"""
import os
import xml.etree.ElementTree as ET
from BioFlow.utils.log_behavior import logger
from BioFlow.utils.IO_Routines import dump_object, undump_object
import BioFlow.main_configs as conf



def zip_dicts(dict1, dict2):
    """
    defines a dictionary update rules that are required for a proper dictionary update

    :param dict1: dictionary #1
    :param dict2: dictionary #2
    """
    for key in dict2.keys():
        if key not in dict1.keys():
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
    """

    # TODO: this needs to be a singleton in order to avoid multiple reloads

    # TODO: more in-depth refactoring:
    #  - it looks like there is a lot of repetitive code because of the similarities in structure
    # and copy-paste style of programming I had back at the time of writing of this module.
    #  - this code is heavily reliant on the string matching and manipulation and case statements
    # are not the best ways of dealing with this.

    # TODO: set the logging.info for all the relevant segments of parsing.

    def __init__(self, path_to_biopax_file=conf.ReactomeBioPax):

        print path_to_biopax_file

        self.tree = ET.parse(path_to_biopax_file)
        self.root = self.tree.getroot()

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

    # all the functions below obey the same pattern:
    #   - recover a tag
    #   - pull the context of a specific tag
    #   - insert it into an appropriated dictionary

    def find_in_root(self, term_name):
        """
        Abstracts the simplest xml foot finding of a path. Currently a test function to be
        integrated in more high-level functions

        :param term_name:
        :return:
        """
        search_pattern = '{http://www.biopax.org/release/biopax-level3.owl#}%s' % term_name
        xml_iterator = self.root.findall(search_pattern)
        print(search_pattern)
        return xml_iterator

    @staticmethod
    def simple_double_tag(sub_dict, _property, term1, term_set2):
        """
        A function that manages dictionary insert

        :param sub_dict: selected_dict[key]
        :param _property:
        :param term1:
        :param term_set2:
        :return:
        """
        if term1 in _property.tag and _property.text:
            for term2 in term_set2:
                if term2 in _property.text:
                    for word in _property.text.split():
                        if term2 in word:
                            sub_dict[term2[:-1]] = word.split(':')[1]

    def parse_tags_with_secondary_switch(self, term_name, tag_to_behavior_dict):
        pass

    def parse_bio_source(self):
        bio_sources_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}BioSource')
        for single_bio_source in bio_sources_xml:
            key = single_bio_source.attrib.values()[0]
            for bio_source_ref_property in single_bio_source:
                if '}name' in bio_source_ref_property.tag:
                    self.BioSources[key] = bio_source_ref_property.text
                    break

    def parse_cellular_locations(self):
        cellular_locations_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}CellularLocationVocabulary')
        for single_cellular_location in cellular_locations_xml:
            key = single_cellular_location.attrib.values()[0]
            for loc_ref in single_cellular_location:
                if '}term' in loc_ref.tag:
                    self.CellularLocations[key] = loc_ref.text

    def parse_sequence_modification_vocabulary(self):
        seq_mod_voc_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}SequenceModificationVocabulary')
        for single_SeqMod in seq_mod_voc_xml:
            key = single_SeqMod.attrib.values()[0]
            for mod_ref in single_SeqMod:
                if '}term' in mod_ref.tag:
                    self.SeqModVoc[key] = mod_ref.text

    def parse_seq_site(self):
        seq_site_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}SequenceSite')
        for single_SeqSite in seq_site_xml:
            key = single_SeqSite.attrib.values()[0]
            for site_ref in single_SeqSite:
                if '}sequencePosition' in site_ref.tag:
                    self.SeqSite[key] = site_ref.text

    # this is the first function that acts differently. After calling the common core, it
    # parses two tags and not one. In addition to that, it calls a more complex function running
    # additional verifications within the tag.

    # maybe a good pattern to refactor that function would be to:
    #  - factor out the common core ({http://www.biopax.org/release/biopax-level3.owl#})
    #  - match specific context in the property tags
    #  form that point our code only uses two variables: the key that was retrieved from the xml
    # tree and
    #  - execute a function that is either a simple lambda to insert key or can be a more elaborate
    #  - static function

    # instead of using an overkill of a static function, we could use a dict ot represent a
    # case-switch with default statement. That would resolve the following patterns:

    #     if 'ENSEMBL:' in dna_ref_property.text:
    #     line = dna_ref_property.text
    #     line = line.split()
    #     for word in line:
    #         if 'ENSEMBL:' in word:
    #             self.DnaRefs[key]['ENSEMBL'] = word.split(':')[1]

    #  if a more complex behavior is encountered, it can be implemented as a separate function and
    #  injected later on.

    def parse_dna_refs(self):
        dna_refs_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}DnaReference')
        for single_Dna_Ref in dna_refs_xml:
            key = single_Dna_Ref.attrib.values()[0]
            self.DnaRefs[key] = {'name': []}
            for dna_ref_property in single_Dna_Ref:
                if '}name' in dna_ref_property.tag and dna_ref_property.text:

                    if 'ENSEMBL:' in dna_ref_property.text:
                        line = dna_ref_property.text
                        line = line.split()
                        for word in line:
                            if 'ENSEMBL:' in word:
                                self.DnaRefs[key]['ENSEMBL'] = word.split(':')[1]
                    else:
                        self.DnaRefs[key]['name'].append(dna_ref_property.text)
                if '}organism' in dna_ref_property.tag and dna_ref_property.attrib.values()[0][
                        1:] != 'BioSource1':
                    self.DnaRefs[key]['organism'] = self.BioSources[
                        dna_ref_property.attrib.values()[0][1:]]

    def parse_rna_refs(self):
        # TODO: refactor, cyclomatic complexity is too high: 14
        rna_refs_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}RnaReference')
        for single_Rna_Ref in rna_refs_xml:
            key = single_Rna_Ref.attrib.values()[0]
            self.RnaRefs[key] = {'name': []}
            for rna_ref_property in single_Rna_Ref:
                if '}name' in rna_ref_property.tag:
                    if 'ENSEMBL:' in rna_ref_property.text:
                        line = rna_ref_property.text
                        line = line.split()
                        for word in line:
                            if 'ENSEMBL:' in word:
                                self.RnaRefs[key]['ENSEMBL'] = word.split(':')[1]
                    else:
                        if 'miRBase:' in rna_ref_property.text:
                            line = rna_ref_property.text
                            line = line.split()
                            for word in line:
                                if 'miRBase:' in word:
                                    self.RnaRefs[key]['miRBase'] = word.split(':')[1]
                        else:
                            if 'EMBL:' in rna_ref_property.text:
                                line = rna_ref_property.text
                                line = line.split()
                                for word in line:
                                    if 'EMBL:' in word:
                                        self.RnaRefs[key]['EMBL'] = word.split(':')[1]
                            else:
                                self.RnaRefs[key]['name'].append(rna_ref_property.text)
                if '}organism' in rna_ref_property.tag and rna_ref_property.attrib.values()[0][
                        1:] != 'BioSource1':
                    self.RnaRefs[key]['organism'] = self.BioSources[
                        rna_ref_property.attrib.values()[0][1:]]

    def parse_small_molecules_refs(self):
        small_molecules_refs_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}SmallMoleculeReference')
        for single_SmallMolecule_Ref in small_molecules_refs_xml:
            key = single_SmallMolecule_Ref.attrib.values()[0]
            self.SmallMoleculeRefs[key] = {'name': []}
            for SmallMolecule_ref_property in single_SmallMolecule_Ref:
                if '}name' in SmallMolecule_ref_property.tag:
                    if 'ChEBI:' in SmallMolecule_ref_property.text:
                        line = SmallMolecule_ref_property.text
                        line = line.split()
                        for word in line:
                            if 'ChEBI:' in word:
                                self.SmallMoleculeRefs[key]['ChEBI'] = word.split(':')[
                                    1].strip(']')
                    else:
                        self.SmallMoleculeRefs[key]['name'].append(
                            SmallMolecule_ref_property.text)
                if '}organism' in SmallMolecule_ref_property.tag \
                        and SmallMolecule_ref_property.attrib.values()[
                        0][1:] != 'BioSource1':
                    self.SmallMoleculeRefs[key]['organism'] = self.BioSources[
                        SmallMolecule_ref_property.attrib.values()[0][1:]]

    def parse_protein_refs(self):
        protein_refs_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}ProteinReference')
        for single_Protein_Ref in protein_refs_xml:
            key = single_Protein_Ref.attrib.values()[0]
            self.ProteinRefs[key] = {'name': []}
            for protein_ref_property in single_Protein_Ref:
                if '}name' in protein_ref_property.tag:
                    if 'UniProt:' in protein_ref_property.text:
                        line = protein_ref_property.text
                        line = line.split()
                        for word in line:
                            if 'UniProt:' in word:
                                self.ProteinRefs[key]['UniProt'] = word.split(':')[1]
                    else:
                        self.ProteinRefs[key]['name'].append(protein_ref_property.text)
                if '}organism' in protein_ref_property.tag and protein_ref_property.attrib.values()[
                        0][1:] != 'BioSource1':
                    self.ProteinRefs[key]['organism'] = self.BioSources[
                        protein_ref_property.attrib.values()[0][1:]]

    def parse_modification_features(self):
        modification_features_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}ModificationFeature')
        for single_ModificationFeature in modification_features_xml:
            key = single_ModificationFeature.attrib.values()[0]
            self.ModificationFeatures[key] = {'ID': key}
            for modification_property in single_ModificationFeature:
                if '}featureLocation' in modification_property.tag:
                    self.ModificationFeatures[key]['location'] = self.SeqSite[
                        modification_property.attrib.values()[0][1:]]
                if '}modificationType' in modification_property.tag:
                    self.ModificationFeatures[key]['modification'] = self.SeqModVoc[
                        modification_property.attrib.values()[0][1:]]

    ##########################################################################
    #
    # Core parses : parsing at the same time objects and their collections
    #
    ##########################################################################
    @staticmethod
    def meta_parser_secondary_loop(local_dict, local_property, is_collection):
        """
        Same approach for all the projects

        :param local_dict: local dictionary where the parser needs to be returned
        :param local_property: if
        :param is_collection: True of we are parsing a collection, false otherwise
        :return:
        """
        if '}cellularLocation' in local_property.tag:
            local_dict['cellularLocation'] = local_property.attrib.values()[0][1:]
        if '}displayName' in local_property.tag:
            local_dict['displayName'] = local_property.text
        if '}name' in local_property.tag:
            local_dict['references']['name'].append(local_property.text)
        if '}memberPhysicalEntity' in local_property.tag:
            is_collection = True
            local_dict['collectionMembers'].append(
                local_property.attrib.values()[0][1:])
        return is_collection

    def parse_dna(self):
        dna_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Dna')
        for single_Dna in dna_xml:
            key = single_Dna.attrib.values()[0]
            local_dict = {'collectionMembers': [], 'references': {'name': []}}
            is_collection = False
            for dna_property in single_Dna:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, dna_property, is_collection)
                if '}entityReference' in dna_property.tag:
                    local_dict['references'] = zip_dicts(
                        self.DnaRefs[
                            dna_property.attrib.values()[0][
                                1:]], local_dict['references'])
            if is_collection:
                self.Dna_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.Dnas[key] = local_dict

    def parse_rna(self):
        rna_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Rna')
        for single_Rna in rna_xml:
            key = single_Rna.attrib.values()[0]
            local_dict = {'collectionMembers': [], 'references': {'name': []}}
            is_collection = False
            for Rna_property in single_Rna:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, Rna_property, is_collection)
                if '}entityReference' in Rna_property.tag:
                    local_dict['references'] = zip_dicts(
                        self.RnaRefs[
                            Rna_property.attrib.values()[0][
                                1:]], local_dict['references'])
            if is_collection:
                self.Rna_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.Rnas[key] = local_dict

    def parse_small_molecules(self):
        small_molecules_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}SmallMolecule')
        for single_SmallMolecule in small_molecules_xml:
            key = single_SmallMolecule.attrib.values()[0]
            local_dict = {'collectionMembers': [], 'references': {'name': []}}
            is_collection = False
            for SmallMolecule_property in single_SmallMolecule:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, SmallMolecule_property, is_collection)
                if '}entityReference' in SmallMolecule_property.tag:
                    local_dict['references'] = zip_dicts(
                        self.SmallMoleculeRefs[
                            SmallMolecule_property.attrib.values()[0][
                                1:]], local_dict['references'])
            if is_collection:
                self.SmallMolecule_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.SmallMolecules[key] = local_dict

    def parse_proteins(self):
        proteins_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Protein')
        for single_Protein in proteins_xml:
            key = single_Protein.attrib.values()[0]
            local_dict = {
                'collectionMembers': [],
                'modification': [],
                'references': {
                    'name': []}}
            is_collection = False
            for Protein_property in single_Protein:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, Protein_property, is_collection)
                if '}entityReference' in Protein_property.tag:
                    local_dict['references'] = zip_dicts(
                        self.ProteinRefs[
                            Protein_property.attrib.values()[0][
                                1:]], local_dict['references'])
                if '}feature' in Protein_property.tag \
                        and 'ModificationFeature' in Protein_property.attrib.values()[
                        0]:
                    local_dict['modification'].append(
                        self.ModificationFeatures[
                            Protein_property.attrib.values()[0][
                                1:]])
            if not local_dict['modification']:
                del local_dict['modification']
            if is_collection:
                self.Protein_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.Proteins[key] = local_dict

    def parse_physical_entities(self):
        physical_entities_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}PhysicalEntity')
        for single_PhysicalEntity in physical_entities_xml:
            key = single_PhysicalEntity.attrib.values()[0]
            local_dict = {'collectionMembers': [], 'references': {'name': []}}
            is_collection = False
            for PhysicalEntity_property in single_PhysicalEntity:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, PhysicalEntity_property, is_collection)
            if is_collection:
                self.PhysicalEntity_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.PhysicalEntities[key] = local_dict

    def parse_complexes(self):
        complex_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Complex')
        for single_Complex in complex_xml:
            key = single_Complex.attrib.values()[0]
            local_dict = {
                'collectionMembers': [],
                'parts': [],
                'references': {
                    'name': []}}
            is_collection = False
            for complex_property in single_Complex:
                is_collection = self.meta_parser_secondary_loop(
                    local_dict, complex_property, is_collection)
                if '}component' in complex_property.tag \
                        and 'Stoichiometry' not in complex_property.tag:
                    local_dict['parts'].append(
                        complex_property.attrib.values()[0][1:])
            if is_collection:
                del local_dict['parts']
                self.Complex_Collections[key] = local_dict
            else:
                del local_dict['collectionMembers']
                self.Complexes[key] = local_dict

    def parse_template_reactions(self):
        template_reaction_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}TemplateReaction')
        for single_TemplateReaction in template_reaction_xml:
            key = single_TemplateReaction.attrib.values()[0]
            local_dict = {'references': {'name': []}}
            for TemplateReaction_property in single_TemplateReaction:
                if '}product' in TemplateReaction_property.tag:
                    local_dict['product'] = TemplateReaction_property.attrib.values()[0][
                        1:]
                if '}displayName' in TemplateReaction_property.tag:
                    local_dict['displayName'] = TemplateReaction_property.text
                if '}name' in TemplateReaction_property.tag:
                    local_dict['references'][
                        'name'] = TemplateReaction_property.text
            self.TemplateReactions[key] = local_dict

    def parse_degradations(self):
        degradation_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Degradation')
        for single_Degradation in degradation_xml:
            key = single_Degradation.attrib.values()[0]
            local_dict = {'references': {'name': []}}
            for Degradation_property in single_Degradation:
                if '}left' in Degradation_property.tag:
                    local_dict['product'] = Degradation_property.attrib.values()[0][
                        1:]
                if '}displayName' in Degradation_property.tag:
                    local_dict['displayName'] = Degradation_property.text
                if '}eCNumber' in Degradation_property.tag:
                    local_dict['references']['eCNumber'] = Degradation_property.text
            self.Degradations[key] = local_dict

    def parse_biochemical_reactions(self):
        biochemical_reaction_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}BiochemicalReaction')
        for single_BiochemicalReaction in biochemical_reaction_xml:
            key = single_BiochemicalReaction.attrib.values()[0]
            local_dict = {'right': [], 'left': [], 'references': {'name': []}}
            for BiochemicalReaction_property in single_BiochemicalReaction:
                if '}left' in BiochemicalReaction_property.tag:
                    local_dict['left'].append(
                        BiochemicalReaction_property.attrib.values()[0][1:])
                if '}right' in BiochemicalReaction_property.tag:
                    local_dict['right'].append(
                        BiochemicalReaction_property.attrib.values()[0][1:])
                if '}displayName' in BiochemicalReaction_property.tag:
                    local_dict['displayName'] = BiochemicalReaction_property.text
                if '}eCNumber' in BiochemicalReaction_property.tag:
                    local_dict['references'][
                        'eCNumber'] = BiochemicalReaction_property.text
            self.BiochemicalReactions[key] = local_dict

    def parse_catalysis(self):

        # TODO: refactor, cyclomatic complexity is too high => just factor out three separate
        # methods. They would be looking almost like clonse of one another.

        catalysis_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Catalysis')
        for single_Catalysis in catalysis_xml:
            key = single_Catalysis.attrib.values()[0]
            self.Catalysises[key] = {}
            for catalysis_property in single_Catalysis:
                if '}controlled' in catalysis_property.tag:
                    self.Catalysises[key]['controlled'] = catalysis_property.attrib.values()[0][
                        1:]
                if '}controller' in catalysis_property.tag:
                    self.Catalysises[key]['controller'] = catalysis_property.attrib.values()[0][
                        1:]
                if '}controlType' in catalysis_property.tag:
                    self.Catalysises[key]['ControlType'] = catalysis_property.text

        template_reaction_regulation_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}TemplateReactionRegulation')
        for single_TRR in template_reaction_regulation_xml:
            key = single_TRR.attrib.values()[0]
            self.Catalysises[key] = {}
            for TRR_property in single_TRR:
                if '}controlled' in TRR_property.tag:
                    self.Catalysises[key]['controlled'] = TRR_property.attrib.values()[0][
                        1:]
                if '}controller' in TRR_property.tag:
                    self.Catalysises[key]['controller'] = TRR_property.attrib.values()[0][
                        1:]
                if '}controlType' in TRR_property.tag:
                    self.Catalysises[key]['ControlType'] = TRR_property.text

        # Both classes above are identical and could be folded ingot a single iterator with
        # iteration concatenation

        # the sub-routine below is different from the two above is that it requires an additiona

        control_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Control')
        for single_Control in control_xml:
            key = single_Control.attrib.values()[0]
            self.Catalysises[key] = {}
            for control_property in single_Control:
                if '}controlled' in control_property.tag \
                        and 'Pathway' not in control_property.attrib.values()[0][1:]:
                    self.Catalysises[key]['controlled'] = control_property.attrib.values()[0][1:]
                if '}controller' in control_property.tag:
                    self.Catalysises[key]['controller'] = control_property.attrib.values()[0][1:]
                if '}controlType' in control_property.tag \
                        and 'BiochemicalReaction' in control_property.attrib.values()[
                        0]:
                    self.Catalysises[key]['ControlType'] = control_property.text

    def parse_modulations(self):
        modulations_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Modulation')
        for single_Modulation in modulations_xml:
            key = single_Modulation.attrib.values()[0]
            self.Modulations[key] = {}
            for modulation_property in single_Modulation:
                if '}displayName' in modulation_property.tag:
                    self.Modulations[key]['displayName'] = modulation_property.text
                if '}controller' in modulation_property.tag:
                    self.Modulations[key]['controller'] = modulation_property.attrib.values()[0][
                        1:]
                if '}controlled' in modulation_property.tag:
                    self.Modulations[key]['controlled'] = self.Catalysises[
                        modulation_property.attrib.values()[0][1:]]['controller']
                if '}controlType' in modulation_property.tag:
                    self.Modulations[key]['controlType'] = modulation_property.text

    def parse_pathways(self):
        pathways_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}Pathway')
        for single_Pathway in pathways_xml:
            key = single_Pathway.attrib.values()[0]
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
                if '}pathwayComponent' in pathway_property.tag \
                        and 'Pathway' in pathway_property.attrib.values()[
                        0]:
                    local_dict['components'].append(
                        pathway_property.attrib.values()[0][1:])
                if '}pathwayOrder' in pathway_property.tag:
                    local_dict['PathwayStep'].append(
                        pathway_property.attrib.values()[0][1:])
            self.Pathways[key] = local_dict

    def parse_pathway_steps(self):
        pathway_steps_xml = self.root.findall(
            '{http://www.biopax.org/release/biopax-level3.owl#}PathwayStep')
        exclude = [
            'Modulation',
            'Control',
            'TemplateReactionRegulation',
            'Catalysis']
        for single_Pathway_step in pathway_steps_xml:
            key = single_Pathway_step.attrib.values()[0]
            local_dict = {'components': [], 'nextStep': []}
            for pathway_property in single_Pathway_step:
                if '}stepProcess ' in pathway_property.tag and not any(
                        x in pathway_property.attrib.values()[0] for x in exclude):
                    local_dict['components'].append(
                        pathway_property.attrib.values()[0][1:])
                if '}nextStep ' in pathway_property.tag:
                    local_dict['nextStep'].append(
                        pathway_property.attrib.values()[0][1:])
            self.PathwaySteps[key] = local_dict

    def parse_all(self):

        self.parse_bio_source()
        self.parse_cellular_locations()
        self.parse_sequence_modification_vocabulary()
        self.parse_seq_site()

        self.parse_dna_refs()
        self.parse_rna_refs()
        self.parse_small_molecules_refs()
        self.parse_protein_refs()
        self.parse_modification_features()

        self.parse_dna()
        self.parse_rna()
        self.parse_small_molecules()
        self.parse_proteins()
        self.parse_physical_entities()
        self.parse_complexes()

        self.parse_template_reactions()
        self.parse_degradations()
        self.parse_biochemical_reactions()
        self.parse_catalysis()

        self.parse_modulations()
        self.parse_pathways()
        self.parse_pathway_steps()

        self.parsed = True

    def get_parse_dicts(self):
        return {
        'BioSources': self.BioSources,
        'CellularLocations': self.CellularLocations,
        'SeqModVoc': self.SeqModVoc,
        'SeqSite': self.SeqSite,
        'Pathways': self.Pathways,
        'PathwaySteps': self.PathwaySteps,

        'DnaRefs': self.DnaRefs,
        'RnaRefs': self.RnaRefs,
        'SmallMoleculeRefs': self.SmallMoleculeRefs,
        'ProteinRefs': self.ProteinRefs,
        'ModificationFeatures': self.ModificationFeatures,

        'Dnas': self.Dnas,
        'Dna_Collections': self.Dna_Collections,
        'Rnas': self.Rnas,
        'Rna_Collections': self.Rna_Collections,
        'SmallMolecules': self.SmallMolecules,
        'SmallMolecule_Collections': self.SmallMolecule_Collections,
        'Proteins': self.Proteins,
        'Protein_Collections': self.Protein_Collections,
        'PhysicalEntities': self.PhysicalEntities,
        'PhysicalEntity_Collections': self.PhysicalEntity_Collections,
        'Complexes': self.Complexes,
        'Complex_Collections': self.Complex_Collections,

        'TemplateReactions': self.TemplateReactions,
        'Degradations': self.Degradations,
        'BiochemicalReactions': self.BiochemicalReactions,
        'Catalysises': self.Catalysises,
        'Modulations': self.Modulations,

        'parsed': self.parsed,
        }

if __name__ == "__main__":
    RP = ReactomeParser()
    RP.parse_all()
    dump_object('../../dumps', RP)
