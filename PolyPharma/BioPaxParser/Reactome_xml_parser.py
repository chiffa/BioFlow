'''
Created on 13 mai 2013

@author: Andrei Kucharavy
'''

import xml.etree.ElementTree as ET
import logging

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

tree = ET.parse('C:\Users\User\Documents\UCSD\Parsing_Reactome\Homo sapiens.owl')
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
        if (subchild.tag,subchild.attrib.keys()[0]) not in Subtypes[CurrentType].keys():
            Subtypes[CurrentType][(subchild.tag,subchild.attrib.keys()[0])]=0
        Subtypes[CurrentType][(subchild.tag,subchild.attrib.keys()[0])]+=1
#        if len(subchild)>0:
#            logging.warning("'further children for '%s','%s'",subchild.tag, subchild.attrib)
#            # Never happens, so this is actually a pretty well formulated rdf

for dattype in Subtypes.keys():
    logging.info("'%s' | '%s'",dattype, TypeCount[dattype])
    for elt in Subtypes[dattype]:
        logging.info(" \t '%s' | '%s' ", elt, Subtypes[dattype][elt])

for EntFeat in root.findall('{http://biopax-level3#}RelationshipTypeVocabulary'):
    logging.info("'%s' | '%s'", EntFeat.tag, EntFeat.attrib)
    for subchild in EntFeat:
        logging.info(" \t '%s' | '%s' | '%s' ", subchild.tag, subchild.attrib, subchild.text)
