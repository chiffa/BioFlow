'''
Created on 13 mai 2013

@author: Andrei Kucharavy
'''

import xml.etree.ElementTree as ET
import logging
from PolyPharma.configs import ReactomeBioPax
import random


# Replace too lengthy IEEE compliand term definitions by something more readable by humans
Reactome_ReadabilityDict={'{http://www.biopax.org/release/biopax-level3.owl#}':''}

def run_xml_Doctor(ReadabilityDict):

    # Updating the readability Dict with what is required for expected functionning if the application

    required_Readability_mapping = {'{http://www.w3.org/1999/02/22-rdf-syntax-ns#}':'rdf:',
                                    '{http://www.w3.org/2002/07/owl#}':'owl:'}

    ReadabilityDict_full = dict(ReadabilityDict.items() + required_Readability_mapping.items())

    # Define logging behavior:

    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='dynamics_full.log',
                        filemode='w')

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    # Parse Reactome BioPax

    tree = ET.parse(ReactomeBioPax)
    root = tree.getroot()

    #TODO: correlation between presences?

    def make_readable(stree):
        """
        Increases readability of a sample output for the final user
        """
        for key in ReadabilityDict_full.keys():
            if key in stree:
                return stree.replace(key,ReadabilityDict_full[key])
        return stree

    def update_Counter_subdict(Counter_subdict, lvl0Dict):
        '''
        Adds a lvl0Dict information to the global counter information
        '''
        Counter_subdict['TotalCount']+=1
        lvl0keys=lvl0Dict.keys()
        lvl0keys.remove('objType')
        for subkey in lvl0keys:
            if subkey not in Counter_subdict.keys():
                Counter_subdict[subkey]={'appearances':0,'count':0}
            Counter_subdict[subkey]['appearances']+=1
            Counter_subdict[subkey]['count']+=lvl0Container[subkey]['count']
            lvl1keys=lvl0Dict[subkey].keys()
            lvl1keys.remove('count')
            for pointer in lvl1keys:
                if pointer not in Counter_subdict[subkey].keys():
                    Counter_subdict[subkey][pointer]={'appearances':0,'count':0}
                Counter_subdict[subkey][pointer]['appearances']+=1
                Counter_subdict[subkey][pointer]['count']+=lvl0Container[subkey][pointer]


    Counter={}

    for child in root:
        lvl0Container={'objType':make_readable(child.tag)}
        for subchild in child:
            lvl1Container=[make_readable(subchild.tag),make_readable(subchild.attrib.keys()[0])]
            if 'resource' in subchild.attrib.keys()[0]: #pointing to a node, a non-terminal node
                # we are no more interested in a type (ressource in the xml from owl), but rather the exact
                # type (Protein, RNA, Small molecule, ...) we have to parse the name to which it points,
                # while getting rid of the identifier number
                referenced=''
                for letter in subchild.attrib.values()[0]:
                    if not letter.isdigit():
                        referenced+=letter
                lvl1Container[1]=referenced
            if not lvl1Container[0] in lvl0Container.keys():
                lvl0Container[lvl1Container[0]]={'count':0}
            lvl0Container[lvl1Container[0]]['count']+=1
            if not lvl1Container[1] in lvl0Container[lvl1Container[0]].keys():
                lvl0Container[lvl1Container[0]][lvl1Container[1]]=0
            lvl0Container[lvl1Container[0]][lvl1Container[1]]+=1
        if 'Ontology' not in lvl0Container['objType']:
            if lvl0Container['objType'] not in Counter.keys():
                Counter[lvl0Container['objType']]={'TotalCount':0}
            update_Counter_subdict(Counter[lvl0Container['objType']],lvl0Container)

    for entry in Counter.keys():
        logging.info("%-30s%s",entry,Counter[entry]['TotalCount'])
        keyslvl1=Counter[entry].keys()
        keyslvl1.remove('TotalCount')
        for key_lvl1 in keyslvl1:
            lvl1_apperance_percentage="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1]['appearances'])/float(Counter[entry]['TotalCount'])*100))
            lvl1_count_per_appearance="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1]['count'])/float(Counter[entry][key_lvl1]['appearances'])))
            logging.info("    %-34s%-10s%s", key_lvl1,lvl1_apperance_percentage+' %',lvl1_count_per_appearance)
            keyslvl2=Counter[entry][key_lvl1].keys()
            keyslvl2.remove('appearances')
            keyslvl2.remove('count')
            if 'rdf:datatype' in keyslvl2:
                keyslvl2.remove('rdf:datatype')
            for key_lvl2 in keyslvl2:
                lvl2_apperance_percentage="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1][key_lvl2]['appearances'])/float(Counter[entry][key_lvl1]['appearances'])*100))
                lvl2_count_per_appearance="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1][key_lvl2]['count'])/float(Counter[entry][key_lvl1][key_lvl2]['appearances'])))
                logging.info("\t%-30s%-8s%s", key_lvl2,'   '+lvl2_apperance_percentage+' %',lvl2_count_per_appearance)
        logging.info("\n")


if __name__ == "__main__":
    run_xml_Doctor(Reactome_ReadabilityDict)