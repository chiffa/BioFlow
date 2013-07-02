'''
Created on Jun 27, 2013

@author: andrei
'''

import requests, json, logging 
import xml.etree.ElementTree as ET


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



def pingPDB(UNIPORT_ACNUM,PDB_ID):
    url='http://www.rcsb.org/pdb/protein/'+UNIPORT_ACNUM
#    payload1={'type':'json'}
    payload2={'type':'json','track':'pdbsites','display':PDB_ID}
#    r1=requests.get(url,params=payload1)
    r2=requests.get(url,params=payload2)
#    if r1.status_code!=200:
#        logging.warning('request 1 with parameter %s failed',UNIPORT_ACNUM)
    if r2.status_code!=200:
        logging.warning('request 2 with parameters %s, %s failed', UNIPORT_ACNUM, PDB_ID)
    return r2.text

def recover_PDB_IDs(UNIPROT_ACNUM):
    url="http://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment"
    payload={'query':'P09211'}
    r=requests.get(url,params=payload)
    root=ET.fromstring(r.text)
    PDB_ID_list=[]
    for child in root:
        for subchild in child:
            if 'dbAccessionId' in subchild.attrib.keys():
                PDB_ID_list.append(subchild.attrib['dbAccessionId'])