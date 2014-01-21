'''
Created on Dec 16, 2013
@author: ank
Contains the Scheme and StrictScheme classes that is used to build neo4j database scheme and 
'''
from __future__ import print_function
from pprint import PrettyPrinter


class SchemeHolder(object):
    '''
    Used to build a neo4j scheme object
    '''

    def __init__(self, params=[]):
        '''
        Constructor
        '''
        self.node_types=[]
        self.relation_types=[]
        self.indexes={} # instance_display_name_<class_>idx: [(class, indexed_option),class_]
        self._ghost_indexes={} # used to access the buffer according to different values # actually a dict of dicts
        self.compile_rules={} # (class_affected_by_the_rule,value_affected_by_the_rule): []
        self.to_add=[]
        self.buffer=[]
        self.to_modify=[]
        self.to_delete=[]
    
    def bind_scheme_dict(self,scheme_dict):
        #TODO: scan all the lvl. dict.variables
        # - look for the options ending in _!: create db-level index for them, used for all the classes inheriting from a given class
        # - look for the options ending in _!!: create class-specific index, where a separate index will be created for each
        # - look for the option ending in _!* or _!!: this means a full-text index will have to be created
        
    
    def instance_completer(self):
        ''' performs the task of computing the attributed defined by the @ symbols'''
        pass
    
    def instance_type_verification(self):
        '''
        In theory used for coercing user-defined attributes to a value types enforced by the scheme.
        currently coerces any element type to a string
        '''
        pass
    
    def local_create(self):
        pass
    
    def local_load(self):
        pass
    
     
0if __name__ == "__main__":
    StS=StrictSchemeGenerator()
    StS.execute()