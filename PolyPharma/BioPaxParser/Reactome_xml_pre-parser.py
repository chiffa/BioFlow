'''
Created on 13 mai 2013

@author: Andrei Kucharavy
'''

import xml.etree.ElementTree as ET
import logging
import random

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='dynamics_full.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

tree = ET.parse('/home/andrei/UCSD/Parsing_Reactome/Homo sapiens.owl')
root = tree.getroot()

#TODO: understand how to perform batch-create with bulbs
#TODO: perform the dictionnaries linking the subtypes of grouping nodes 
# ATTENTION: if there is a feature, it should be reagarded as a 

Subtypes={}
TypeCount={}

for child in root:
    CurrentType=child.tag
    if CurrentType not in Subtypes.keys():
        Subtypes[CurrentType]={}
        TypeCount[CurrentType]=0
    TypeCount[CurrentType]+=1
    for subchild in child:
        refference=""
        if 'resource' in subchild.attrib.keys()[0]:
            for letter in subchild.attrib.values()[0]:
                if not letter.isdigit():
                    refference+=letter
        else : refference=subchild.attrib.keys()[0]
        if (subchild.tag,refference) not in Subtypes[CurrentType].keys():
            Subtypes[CurrentType][(subchild.tag,refference)]=0
        Subtypes[CurrentType][(subchild.tag,refference)]+=1
#        if len(subchild)>0:
#            logging.warning("'further children for '%s','%s'",subchild.tag, subchild.attrib)
#            # Never happens, so this is actually a pretty well formulated rdf


for dattype in Subtypes.keys():
    logging.info("'%s' | '%s'",dattype, TypeCount[dattype])
    for elt in Subtypes[dattype]:
        logging.info(" \t '%s' | '%s' ", elt, "{0:.2f}".format(float(Subtypes[dattype][elt])/float(TypeCount[dattype])*100.0))
    logging.info('\n')
    Sample=root.findall(dattype)
    random.shuffle(Sample)
    for EntFeat in Sample[:5]:
        logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
        for subchild in EntFeat:
            logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
    logging.info('<=======================================================>\n\n\n')
    


Sample=root.findall('{http://biopax-level3#}BioSource')
for EntFeat in Sample:
    logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
    for subchild in EntFeat:
        logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
logging.info('<=======================================================>\n\n\n')