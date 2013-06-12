'''
Created on 13 mai 2013

@author: Andrei Kucharavy
'''

import xml.etree.ElementTree as ET
import logging
import random

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

tree = ET.parse('/home/andrei/UCSD/Parsing_Reactome/Homo sapiens.owl')
root = tree.getroot()

#TODO: understand how to perform batch-create with bulbs
#TODO: perform the dictionnaries linking the subtypes of grouping nodes 
# ATTENTION: if there is a feature, it should be reagarded as a 


def make_readable(stree):
    txt1='{http://www.biopax.org/release/biopax-level3.owl#}'
    if txt1 in stree:
        return stree.replace(txt1,'')
    txt2='{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    if txt2 in stree:
        return stree.replace(txt2,'rdf:')
    txt3='{http://www.w3.org/2002/07/owl#}'
    if txt3 in stree:
        return stree.replace(txt3,'owl:')
    return stree

def update_Counter_subdict(Counter_subdict,lvl0Dict):
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
            # we are no more interested in a type (ressource in the xml from owl), but rather the exact type (Protein, RNA, Small molecule, ...)
            # we have to parse the name to which it points, while getting rid of the identifier number
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
    logging.info("%-40s%s",entry,Counter[entry]['TotalCount'])
    keyslvl1=Counter[entry].keys()
    keyslvl1.remove('TotalCount')
    for key_lvl1 in keyslvl1:
        lvl1_apperance_percentage="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1]['appearances'])/float(Counter[entry]['TotalCount'])*100))
        lvl1_count_per_appearance="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1]['count'])/float(Counter[entry][key_lvl1]['appearances'])))
        logging.info("\t  %-40s%-10s%s", key_lvl1,lvl1_apperance_percentage+' %',lvl1_count_per_appearance)
        keyslvl2=Counter[entry][key_lvl1].keys()
        keyslvl2.remove('appearances')
        keyslvl2.remove('count')
        if 'rdf:datatype' in keyslvl2:
            keyslvl2.remove('rdf:datatype')
        for key_lvl2 in keyslvl2:
            lvl2_apperance_percentage="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1][key_lvl2]['appearances'])/float(Counter[entry][key_lvl1]['appearances'])*100))
            lvl2_count_per_appearance="{:>10}".format("{0:.2f}".format(float(Counter[entry][key_lvl1][key_lvl2]['count'])/float(Counter[entry][key_lvl1][key_lvl2]['appearances'])))
            logging.info("\t\t    %-30s%-9s%s", key_lvl2,lvl2_apperance_percentage+' %',lvl2_count_per_appearance)
    logging.info("\n")

#TODO: correlation between presences?
#'''
#<=========================================================================================================================================>
#'''


# Subtypes={}
# Subtypes2={}
# TypeCount={}
# 
# for child in root:
#     CurrentType=child.tag
#     if CurrentType not in Subtypes.keys():
#         Subtypes[CurrentType]={}
#         Subtypes2[CurrentType]={}
#         TypeCount[CurrentType]=0
#     TypeCount[CurrentType]+=1
#     RunningDict=set()
#     for subchild in child:
#         refference=""
#         if 'resource' in subchild.attrib.keys()[0]:
#             for letter in subchild.attrib.values()[0]:
#                 if not letter.isdigit():
#                     refference+=letter
#         else : refference=subchild.attrib.keys()[0]
#         pair=(subchild.tag,refference)
#         if pair not in Subtypes[CurrentType].keys():
#             Subtypes[CurrentType][pair]=0
#             Subtypes2[CurrentType][pair]=0
#         Subtypes[CurrentType][pair]+=1
#         RunningDict.add(pair)
#     for pair in RunningDict:
#         Subtypes2[CurrentType][pair]+=1
#         
#     # Wait, this is logical that the coun2 would be superior, since the signature is encountered several times in different contexts. Like left: reagent1(protein) and left: reagent2(RNA) are both pretty valid
#     # Add a global add, so that only one iteration is performed each time there is a subchild present
#     
# #        if len(subchild)>0:
# #            logging.warning("'further children for '%s','%s'",subchild.tag, subchild.attrib)
# #            # Never happens, so this is actually a pretty well formulated rdf
# 
# Restructured={}
# 
# for objtype in Subtypes.keys(): #they are unique, so Subtypes2.keys()==Subtypes.keys()=True
#     Restructured[objtype]={}
#     for daattype in Subtypes[objtype]:
#         if not daattype[0] in Restructured[objtype].keys():
#             Restructured[objtype][daattype[0]]={"count1":0,"count2":0}
#         Restructured[objtype][daattype[0]]["count1"]+=Subtypes[objtype][daattype]
#         logging.debug("debug: added to count1 %s based on %s, %s", Subtypes[objtype][daattype],objtype,daattype)
#         Restructured[objtype][daattype[0]]["count2"]+=Subtypes2[objtype][daattype]
#         logging.debug("debug: added to count2 %s based on %s, %s", Subtypes2[objtype][daattype],objtype,daattype)
#         Restructured[objtype][daattype[0]][daattype[1]]=(Subtypes[objtype][daattype], Subtypes2[objtype][daattype])
#         logging.debug("\t<<")
#     logging.debug("\n")
# 
# for objtype in Restructured.keys():
#     logging.info("'%s' | '%s'", objtype, TypeCount[objtype])
#     for dattype in Restructured[objtype].keys():
#         logging.info(" \t '%s' | present in '%s' percent | on average: '%s' items", dattype, "{0:.2f}".format(float(max(Restructured[objtype][dattype]["count2"],1))/float(TypeCount[objtype])*100.0), "{0:.2f}".format(float(Restructured[objtype][dattype]["count1"])/float(max(Restructured[objtype][dattype]["count2"],1)))) 
#         for subdattype in Restructured[objtype][dattype].keys():
#             if subdattype!="count1" and subdattype!="count2":
#                 logging.info(" \t\t '%s' | present in '%s' percent | on average: '%s' items", subdattype, "{0:.2f}".format(float(max(Restructured[objtype][dattype][subdattype][1],1))/float(max(Restructured[objtype][dattype]["count2"],1))*100.0), "{0:.2f}".format(float(Restructured[objtype][dattype][subdattype][0])/float(max(Restructured[objtype][dattype][subdattype][1],1)))) 
#     logging.info('<=======================================================>\n\n')
# 
# 
# for dattype in Subtypes.keys():
#     logging.info("'%s' | '%s'",dattype, TypeCount[dattype])
#     for elt in Subtypes[dattype]:
#         logging.info(" \t '%s' | present in '%s'% | on average: '%s' items", elt, "{0:.2f}".format(float(max(Subtypes2[dattype][elt],1))/float(TypeCount[dattype])*100.0), "{0:.2f}".format(float(Subtypes[dattype][elt])/float(max(Subtypes2[dattype][elt],1)))) 
# #     logging.info('\n')
# #     Sample=root.findall(dattype)
# #     random.shuffle(Sample)
# #     for EntFeat in Sample[:5]:
# #         logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
# #         for subchild in EntFeat: 
# #             logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
#     logging.info('<=======================================================>\n\n\n')
#     
# 
# 
# Sample=root.findall('{http://biopax-level3#}BioSource')
# for EntFeat in Sample:
#     logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
#     for subchild in EntFeat:
#         logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
# logging.info('<=======================================================>\n\n\n')