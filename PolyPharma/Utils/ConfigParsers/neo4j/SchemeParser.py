'''
Created on Dec 16, 2013
@author: ank
A class specifically designed for communication with the 
'''

from __future__ import print_function
import ConfigParser
import os
from pprint import PrettyPrinter
from itertools import compress
import warnings

rootdir=os.path.abspath(os.path.join(os.path.dirname(__file__),'../schemes/'))
shortnames=['MainScheme','NoteScheme']
schemefiles=[rootdir+'/'+name+'.ini' for name in shortnames ]
ErrorString="<broken inheritence!>"
baseclasses=['node', 'adjacency', 'relation']
# Additional levels of restriction, we'll leave them out for now
allowed_container_types=["_String","_Float","_Int"]
allowed_class_types=["_Index","_Node","_Relation"]

def generate_scheme():
    '''
    '''
    scheme = ConfigParser.SafeConfigParser()

    scheme.set('DEFAULT', "*_inherits_from", ErrorString)       # Is used to infer type. In case of classType is default, inferred classtype is used, otherwise a warning is issued
    scheme.set('DEFAULT', "Instance_Display_Name_!", "@self.class_name+@self._instance_ID")
    scheme.set('DEFAULT', "*_instance_ID_!", "*_String")
    scheme.add_section('NODE')
    scheme.set('NODE', "*_class_type", 'Node')
    scheme.set('NODE', "*_inherits_from", '@BaseClass')  
    scheme.add_section('RELATION')
    scheme.set('RELATION', "*_class_type", 'Relation')
    scheme.set('RELATION', "*_inherits_from", '@BaseClass')      # Is used to infer type. In case of classType is default, inferred classtype is used, otherwise a warning is issued    
    scheme.add_section('ADJACENCY')                              #Will only be used for the strict schemes; RESPECTS INHERITENCE MODELS
    scheme.set('ADJACENCY', "*_class_type", 'Adjacency')
    scheme.set('ADJACENCY', "*_inherits_from", '@BaseClass')
    scheme.set('ADJACENCY', "*_permitted_connects", '@node, @relation')
    with open(schemefiles[0],'w') as configfile:
        scheme.write(configfile)

def parse_scheme(path_to_scheme,schemeStorageObject=""):
    ''' '''

    def improved_read(path):
        '''
        Will throw an IOError in case it can't read a file
        '''
        scheme_ini=ConfigParser.SafeConfigParser()
        rfs=scheme_ini.read(path)
        if rfs==[]:
            raise IOError('cannot load '+path)
        MainDict={}
        for section in scheme_ini.sections():
            section_name="_".join([elt.lower() for elt in section.split('_')])
            MainDict[section_name]={}
            for option in scheme_ini.options(section):
                MainDict[section_name][option]=scheme_ini.get(section, option).lower()
        return MainDict

    def invert_inheritence(arg_Inheritant2Inherited):
        '''inverts the inheritence dictionary'''
        Inerited2Inheritant={}
        for key, value in arg_Inheritant2Inherited.iteritems():
            if value not in Inerited2Inheritant.keys():
                Inerited2Inheritant[value]=[]
            Inerited2Inheritant[value].append(key)
        return Inerited2Inheritant
    
    def check_min_config(scheme_dict):        
        
        def strict_lvl2_flatten(dico):
            '''
            Flattens dicts up to level 2 deep, assuming all the values in the global dictionary are nested dictionaries
            and all the keys and values are strings
            '''
            FlatDict={}
            for superkey in dico.keys():
                for subkey in dico[superkey].keys():
                    FlatDict[superkey+"/"+subkey]=dico[superkey][subkey]
            return FlatDict
        
        def ref_scheme_error_print(erred_Dict_List,flatenend_actual_dict):
            ''' prints beautifully the error pairs (as much as an error can be beautiful) '''
            returnString="\n\n"
            for elt in erred_Dict_List:
                returnString+=elt[0].split('/')[0].upper()+ " = "+elt[1]+ " expected\n"
                coelt=flatenend_actual_dict.get(elt[0])
                try:
                    returnString+=elt[0].split('/')[0].upper()+" = "+coelt+" found\n\n"
                except TypeError:
                    returnString.append("nothing found\n")
            return returnString
        
        ref_scheme=improved_read(schemefiles[0])
        d1,d2=(strict_lvl2_flatten(ref_scheme),strict_lvl2_flatten(scheme_dict))
        if not all(item in scheme_dict.items() for item in ref_scheme.items()):
            errlist=list(compress(d1.items(),[item not in d2.items() for item in d1.items()]))
            errstr="Reference classes have been altered. Please return Default, Node, Relation and Adjacency classes in their initial states:"+ref_scheme_error_print(errlist,d2)
            raise Exception(errstr)
        if ErrorString in d2.values():
            errstr="A class that inherits from no-one has been found. He feels lonely. Please correct this for: "+inheritance_error_print(d2)
            raise Exception(errstr)
        return d2

    def check_connexity(flattened_dict):
        ''' in addition to checking that no loops are present, also checks that none of the classes is attached to the baseclass'''

        def inheritance_error_print(flattened_actual_dict):
            ''' prints all the nodes that cause inheritence error'''
            returnString="\n"
            for composite_key, value in flattened_actual_dict.iteritems():
                if value==ErrorString:
                    returnString+=composite_key+"\n"
            return returnString
        
        def reduce_inheritance_dictionnary(arg_Inh_Dict):
            Inh_Dict=arg_Inh_Dict.copy()
            while True:
                prevDict=Inh_Dict.copy()
                leaf_classes=set(Inh_Dict.keys()).difference(set(Inh_Dict.values()))
                for cls in leaf_classes:
                    if Inh_Dict[cls]!="baseclass":
                        del Inh_Dict[cls]
                if prevDict==Inh_Dict:
                    break
            return Inh_Dict
    
        def build_inheritance(flattened_actual_dict):
            ''' '''
            Inheritant2Inherited={}
            for key, value in flattened_actual_dict.iteritems():
                if "*_inherits_from" in key:
                    Inheritant2Inherited[key.split('/')[0].lower()]=value.strip('@')
            for cls in baseclasses:
                del Inheritant2Inherited[cls]
            return Inheritant2Inherited 
        
        InhDict=build_inheritance(flattened_dict)
        Reduced_Class_Dict=reduce_inheritance_dictionnary(InhDict)
        if len(Reduced_Class_Dict)==0:
            return InhDict
        Inerited2Inheritant=invert_inheritence(Reduced_Class_Dict)
        if "baseclass" in Inerited2Inheritant.keys():
            errmsg="Non-base classes inherit from the baseclass: \n"+Inerited2Inheritant["baseclass"]
            raise Exception(errmsg)
        else: 
            errmsg="Circular inheritance for the following objects: \n"+Inerited2Inheritant.keys().__str__()
            raise Exception(errmsg)
        #TODO:feature, add verification of different inheritance loops  
   
    def inflate_classes(scheme_dict,Inheritant2Inherited):
        
        def copy_with_proper_overrides(arg_raw_class,arg_inflated_father_class):
            # fine-grained override section
            raw_class=arg_raw_class.copy()
            inflated_father_class=arg_inflated_father_class.copy()
            raw_keyset=set(raw_class.keys())
            father_keyset=set(inflated_father_class.keys())
            interset=raw_keyset.intersection(father_keyset)
            for key in interset:
                if key=='*_class_type' and inflated_father_class[key]!=raw_class[key]:
                    override_msg='Overriding a different type declaration:\n '
                    override_msg+=arg_raw_class['*_class_name']+' is declared as '+raw_class[key]
                    override_msg+=',\n but inherited from '+arg_inflated_father_class['*_class_name']+' which is of class '+inflated_father_class[key]
                    warnings.warn(override_msg)
                    del raw_class[key]
                else:
                    del inflated_father_class[key]
            inflated_father_class.update(raw_class)
            return inflated_father_class
            
        def recursive_cls_inflate(inherited_cls_dict,Reverse_Inh_Dict,collection_dict,raw_scheme):
            ''' inflates recursively'''
            # Stop recursion
            if inherited_cls_dict not in Reverse_Inh_Dict.keys():
                return collection_dict
            else:
                for cls in Reverse_Inh_Dict[inherited_cls_dict]:
                    raw_scheme[cls]['*_class_name']=cls
                    del raw_scheme[cls]['*_inherits_from']
                    collection_dict[cls]=copy_with_proper_overrides(raw_scheme[cls], collection_dict[inherited_cls_dict])
                    collection_dict=recursive_cls_inflate(cls, Reverse_Inh_Dict, collection_dict, raw_scheme)
                return collection_dict
        
        ReverseInheritence=invert_inheritence(Inheritant2Inherited)
        collection_dict={}
        for cls in baseclasses:
            scheme_dict[cls]['*_class_name']=cls
            del scheme_dict[cls]['*_inherits_from']
            collection_dict[cls]=scheme_dict[cls]
            recursive_cls_inflate(cls, ReverseInheritence,  collection_dict, scheme_dict)
        
        return collection_dict, ReverseInheritence

    def process():
        scheme_dict=improved_read(path_to_scheme)
        flat_dict = check_min_config(scheme_dict)
        Inheritant2Inherited=check_connexity(flat_dict)
        return inflate_classes(scheme_dict,Inheritant2Inherited)

    # TODO: Edit to actually get the inheritence tree for adjacency-restricted operations
    return process()[0]

if __name__ == "__main__":
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)
    generate_scheme()
    pp=PrettyPrinter(indent=4)
    pp.pprint(parse_scheme(schemefiles[1]))