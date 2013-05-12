'''
Created on Mar 27, 2013

@author: akucahravy
'''

# TODO: in the protein and GO reporting modules filter out the GOs and proteins 
# without information circulation, eventually throw prots as an additional output

class Protein():
    '''
    Wraps the protein object and the relation it has with GOs and other nodes
    '''
    Flow_With_All=0.0 #Flow with all the unmasked proteins
    Importance=1.0
    GOList=[] # List of all the GOs annotating the Protein
    Flow_to_Prots={} #Dict of type {ProtName : value} representing the total InfFlow
                        # Between two proteins
    Tabulist=[] 
    
    


class MyClass(object):
    '''
    classdocs
    '''


    def __init__(selfparams):
        '''
        Constructor
        '''
        