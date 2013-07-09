'''
Created on Jul 5, 2013

@author: andrei
'''

from Utils.GO_Structure_Parser import GO_Terms, GO_Terms_Structure
from Utils.UNIPROT_Parser import Uniprot
from neo4j_Declarations.Graph_Declarator import DatabaseGraph

GODict={} # Stores relations between GO IDs and the objects in the neo4j database
UniprotDict={} # Stores relations between the SWISSPROT UNIPROT IDs and the neo4j database objects

def import_GOs():
    # Create Terms
    for GO_Term in GO_Terms.keys():
        primary=DatabaseGraph.GOTerm.create(ID=GO_Terms[GO_Term]['id'],Name=GO_Terms[GO_Term]['name'],displayName=GO_Terms[GO_Term]['name'],Namespace=GO_Terms[GO_Term]['namespace'],Definition=GO_Terms[GO_Term]['def'])
        GODict[GO_Term]=primary
    # Create the structure between them:
    for relkey in GO_Terms_Structure.keys():
        primary=GODict[relkey]
        secondary=GODict[GO_Terms_Structure[relkey][1]]
        if GO_Terms_Structure[relkey][0]=='is_a':
            DatabaseGraph.is_a_go.create(primary,secondary)
        if GO_Terms_Structure[relkey][0]=='part_of':
            DatabaseGraph.is_part_of_go.create(primary,secondary)
        if 'regul' in GO_Terms_Structure[relkey][0]:
            if GO_Terms_Structure[relkey][0]=='positively_regulates':
                DatabaseGraph.is_regulant.create(primary,secondary,controlType='ACTIVATES',ID=str('GO'+primary.ID+secondary.ID))
            if GO_Terms_Structure[relkey][0]=='negatively_regulates':
                DatabaseGraph.is_regulant.create(primary,secondary,controlType='INHIBITS',ID=str('GO'+primary.ID+secondary.ID))
            else:
                DatabaseGraph.is_regulant.create(primary,secondary,ID=str('GO'+primary.ID+secondary.ID))

def import_UNIPROTS():
    for CH_PROT_ID in Uniprot.keys():
        #Create uniprot terms
        primary=DatabaseGraph.UNIPORT.create(ID=CH_PROT_ID,displayName=Uniprot[CH_PROT_ID]['Names']['Full'])
        UniprotDict[CH_PROT_ID]=primary
        # Insert references to GOs
        for GO_Term in Uniprot[CH_PROT_ID]['GO']:
            secondary=GODict[GO_Term]
            DatabaseGraph.is_go_annotation.create(primary,secondary)
        # Find intersection between logged acnums and acnums pointed out by the Reactome.org. Create a direct bridge between a reactome ID and a SWISSPROT ID
        for acnum in Uniprot[CH_PROT_ID]['Acnum']:
            generator=DatabaseGraph.AnnotNode.index.lookup(payload=acnum)
            lst=[]
            for elt in generator:
                if elt!=None:
                    ID=str(elt).split('/')[-1][:-1]
                    lst.append(ID)
            if len(lst)>0:
                print lst
                for currentID in lst:
                    CurrentNode=DatabaseGraph.vertices.get(currentID)
                    seekedproteinGen=CurrentNode.bothV()
                    lst2=[]
                    for prot in seekedproteinGen:
                        if prot!=None:
                            ID=str(prot).split('/')[-1][:-1]
                            protObj=DatabaseGraph.vertices.get(ID)
                            lst2.append(protObj)
                    for elt in lst2:
                        DatabaseGraph.is_same(primary,elt)
                            
    
    
    
    


        