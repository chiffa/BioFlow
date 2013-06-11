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

Subtypes={}
Subtypes2={}
TypeCount={}

for child in root:
    CurrentType=child.tag
    if CurrentType not in Subtypes.keys():
        Subtypes[CurrentType]={}
        Subtypes2[CurrentType]={}
        TypeCount[CurrentType]=0
    TypeCount[CurrentType]+=1
    RunningDict=set()
    for subchild in child:
        refference=""
        if 'resource' in subchild.attrib.keys()[0]:
            for letter in subchild.attrib.values()[0]:
                if not letter.isdigit():
                    refference+=letter
        else : refference=subchild.attrib.keys()[0]
        pair=(subchild.tag,refference)
        if pair not in Subtypes[CurrentType].keys():
            Subtypes[CurrentType][pair]=0
            Subtypes2[CurrentType][pair]=0
        Subtypes[CurrentType][pair]+=1
        RunningDict.add(pair)
    for pair in RunningDict:
        Subtypes2[CurrentType][pair]+=1
        
    # Add a global add, so that only one iteration is performed each time there is a subchild present
    
#        if len(subchild)>0:
#            logging.warning("'further children for '%s','%s'",subchild.tag, subchild.attrib)
#            # Never happens, so this is actually a pretty well formulated rdf

Restructured={}

for objtype in Subtypes.keys():
    Restructured[objtype]={}
    for daattype in Subtypes[objtype]:
        if not daattype[0] in Restructured[objtype].keys():
            Restructured[objtype][daattype[0]]={"count1":0,"count2":0}
        Restructured[objtype][daattype[0]]["count1"]+=Subtypes[objtype][daattype]
        Restructured[objtype][daattype[0]]["count2"]+=Subtypes2[objtype][daattype]
        Restructured[objtype][daattype[0]][daattype[1]]=(Subtypes[objtype][daattype], Subtypes2[objtype][daattype])

for objtype in Restructured.keys():
    logging.info("'%s' | '%s'", objtype, TypeCount[objtype])
    for dattype in Restructured[objtype].keys():
        logging.info(" \t '%s' | present in '%s' percent | on average: '%s' items", dattype, "{0:.2f}".format(float(max(Restructured[objtype][dattype]["count2"],1))/float(TypeCount[objtype])*100.0), "{0:.2f}".format(float(Restructured[objtype][dattype]["count1"])/float(max(Restructured[objtype][dattype]["count2"],1)))) 
        for subdattype in Restructured[objtype][dattype].keys():
            if subdattype!="count1" and subdattype!="count2":
                logging.info(" \t\t '%s' | present in '%s' percent | on average: '%s' items", subdattype, "{0:.2f}".format(float(max(Restructured[objtype][dattype][subdattype][1],1))/float(max(Restructured[objtype][dattype]["count2"],1))*100.0), "{0:.2f}".format(float(Restructured[objtype][dattype][subdattype][0])/float(max(Restructured[objtype][dattype][subdattype][1],1)))) 
    logging.info('<=======================================================>\n\n')


for dattype in Subtypes.keys():
    logging.info("'%s' | '%s'",dattype, TypeCount[dattype])
    for elt in Subtypes[dattype]:
        logging.info(" \t '%s' | present in '%s'% | on average: '%s' items", elt, "{0:.2f}".format(float(max(Subtypes2[dattype][elt],1))/float(TypeCount[dattype])*100.0), "{0:.2f}".format(float(Subtypes[dattype][elt])/float(max(Subtypes2[dattype][elt],1)))) 
#     logging.info('\n')
#     Sample=root.findall(dattype)
#     random.shuffle(Sample)
#     for EntFeat in Sample[:5]:
#         logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
#         for subchild in EntFeat: 
#             logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
    logging.info('<=======================================================>\n\n\n')
    


Sample=root.findall('{http://biopax-level3#}BioSource')
for EntFeat in Sample:
    logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
    for subchild in EntFeat:
        logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
logging.info('<=======================================================>\n\n\n')