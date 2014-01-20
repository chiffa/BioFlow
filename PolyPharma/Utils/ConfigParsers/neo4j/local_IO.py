'''
Created on Dec 15, 2013
@author: ank
Deals with the local IO for the neo4j database
'''

import src.ConfigParser as cfgPrs
import logging
import os

##########  Configures the options for logging and JVM parameters  ############
cfg1, cfg2 = cfgPrs.parse_configs()  
for key in cfg2['JVM']:
    os.environ[key.upper()]=cfg2['JVM'][key]
dico = {'DEBUG':logging.DEBUG, 'INFO':logging.INFO,'WARNING':logging.WARNING,
         'ERROR':logging.ERROR, 'CRITICAL':logging.CRITICAL}
try:
    lvl=dico[cfg2['LOGGER']['level']]
except KeyError:
    lvl=logging.WARNING
logging.basicConfig(filename=cfg2['LOGGER']['output_file'],level=lvl)
###############################################################################

from neo4j import GraphDatabase

class schemeObject():
    '''
    Represents a scheme for communication with a neo4j instance
    Loads from a configFile
    '''

class Buffer_Container(object):
    ''' Buffering class to avoid unneccesary IO with the noe4j database'''
    
    def __init__(self,schemeObject):
        self.dict={}
        pass
    

class Connector(object):
    ''' Connects and configures a local neo4j instance and a Container object'''
    
    def __init__(self, local_server_path=cfg1['TEST']["local_neo4j"],schemeObject):
        '''
        neo is a true neo4j database.
        '''
        self.neo=GraphDatabase(local_server_path)
        self.Buffer=Buffer_Container(schemeObject)

    def load_from_scheme(self,schemeObject):
        pass
        
    def create_from_scheme(self,schemeObject):
        pass
             
      

#Containers contain indexes and objects that have been loaded from the scheme-defined objects
 
# def create_scheme(text_file):
#     '''
#     creates a scheme file from a text file containing 
#     '''
# 
# def load_from_scheme(neo, schemeObject, container):
#     with neo.transaction:
#         prot_idx=neo.node.indexes.get('proteins')
#         gos_idx=neo.node.indexes.get('gos')
#         for rel in neo.reference_node.GOS:
#             gos = rel.end
#         for rel in neo.reference_node.PROTEINS:
#             proteins = rel.end
# 
# def create_from_scheme(neo, schemeObject,container):
#     with neo.transaction:
#         proteins=neo.node()
#         gos=neo.node()
#         neo.reference_node.PROTEINS(proteins)
#         neo.reference_node.GOS(gos)
#         prot_idx=neo.node.indexes.create('proteins')
#         gos_idx=neo.node.indexes.create('gos')