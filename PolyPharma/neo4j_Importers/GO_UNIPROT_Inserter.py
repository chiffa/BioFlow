'''
Created on Jul 5, 2013

@author: andrei
'''

from Utils.GO_Structure_Parser import GO_Terms, GO_Terms_Structure
from Utils.UNIPROT_Parser import Uniprot
from neo4j_Declarations.Graph_Declarator import DatabaseGraph
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='Go_UP_insert_log.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

GODict={} # Stores relations between GO IDs and the objects in the neo4j database
UniprotDict={} # Stores relations between the SWISSPROT UNIPROT IDs and the neo4j database objects

def import_GOs():
    # Create Terms
    leng=len(GO_Terms.keys())
    i=0
    for GO_Term in GO_Terms.keys():
        i+=1
        logging.debug('GO %s', str("{0:.2f}".format(float(i)/float(leng)*100)))
        primary=DatabaseGraph.GOTerm.create(ID=GO_Terms[GO_Term]['id'],Name=GO_Terms[GO_Term]['name'],displayName=GO_Terms[GO_Term]['name'],Namespace=GO_Terms[GO_Term]['namespace'],Definition=GO_Terms[GO_Term]['def'])
        GODict[GO_Term]=primary
    # Create the structure between them:
    leng=len(GO_Terms_Structure)
    i=0
    for relation in GO_Terms_Structure:
        i+=1
        logging.debug('rel %s', str("{0:.2f}".format(float(i)/float(leng)*100)))
        primary=GODict[relation[0]]
        secondary=GODict[relation[2]]
        Type=relation[1]
        if Type=='is_a':
            DatabaseGraph.is_a_go.create(primary,secondary)
        if Type=='part_of':
            DatabaseGraph.is_part_of_go.create(primary,secondary)
        if 'regul' in Type:
            if Type=='positively_regulates':
                DatabaseGraph.is_regulant.create(primary,secondary,controlType='ACTIVATES',ID=str('GO'+primary.ID+secondary.ID))
            if Type=='negatively_regulates':
                DatabaseGraph.is_regulant.create(primary,secondary,controlType='INHIBITS',ID=str('GO'+primary.ID+secondary.ID))
            else:
                DatabaseGraph.is_regulant.create(primary,secondary,ID=str('GO'+primary.ID+secondary.ID))

def loadAcnumDecs():
    import Reactome_org_parser as RP
    RP_pID={}
    for proteinID in RP.Proteins.keys():
        if 'UniProt' in RP.Proteins[proteinID]['references'].keys():
            RP_pID[RP.Proteins[proteinID]['references']['UniProt']]=proteinID
    Acnu2ProtObj={}
    for key in RP_pID.keys():
        Acnu2ProtObj[key]=[]
        generator=DatabaseGraph.Protein.index.lookup(ID=RP_pID[key])
        if generator!=None:
            for elt in generator:
                ID=str(elt).split('/')[-1][:-1]
                ProtObj=DatabaseGraph.vertices.get(ID)
                Acnu2ProtObj[key].append(ProtObj)
    logging.debug('2 %s', len(Acnu2ProtObj))
    return Acnu2ProtObj


def getExistingAcnums():
    # Try to retrieve acnums from the neo4j database, if fails uses the parsed relations
    generator=DatabaseGraph.AnnotNode.index.lookup(ptype='UniProt')
    if generator==None:
        return loadAcnumDecs()
    AcnumList={} #acnum to AnnotNode
    for elt in generator:
        if elt!=None:
            ID=str(elt).split('/')[-1][:-1]
            AnnotObj=DatabaseGraph.vertices.get(ID)
            AcnumList[str(AnnotObj.payload)]=AnnotObj
    if len(AcnumList)<10:
        return loadAcnumDecs()
    ReactProtList={}
    for key in AcnumList.keys():
        ReactProtGen=AcnumList[key].bothV()
        ReactProtList[key]=[]
        if ReactProtGen!=None:
            for vertex in ReactProtGen:
                if vertex!=None:
                    ReactProtList[key].append(vertex)
    logging.debug('1 %s', len(ReactProtList))
    return ReactProtList
    

def manage_acnums(acnum,Acnum2RProts):
    if acnum not in Acnum2RProts.keys():
        return []
    if acnum in Acnum2RProts.keys():
        return Acnum2RProts[acnum]

def import_UNIPROTS():
    Acnums2RProts=getExistingAcnums()
    i=0
    j=0
    leng=len(Acnums2RProts.keys())
    for CH_PROT_ID in Uniprot.keys():
        set1=set(Uniprot[CH_PROT_ID]['Acnum'])
        set2=set(Acnums2RProts.keys())
        if not set1.isdisjoint(set2):
            i+=1
            logging.debug('UNIPROT %s', str("{0:.2f}".format(float(i)/float(leng)*100)))
            #Create uniprot terms
            primary=DatabaseGraph.UNIPORT.create(ID=CH_PROT_ID,displayName=Uniprot[CH_PROT_ID]['Names']['Full'])
            UniprotDict[CH_PROT_ID]=primary
            # Insert references to GOs
            for GO_Term in Uniprot[CH_PROT_ID]['GO']:
                if GO_Term in GODict.keys():
                    secondary=GODict[GO_Term]
                    DatabaseGraph.is_go_annotation.create(primary,secondary)
            # Find intersection between logged acnums and acnums pointed out by the Reactome.org. Create a direct bridge between a reactome ID and a SWISSPROT ID
            for acnum in Uniprot[CH_PROT_ID]['Acnum']:
                proteins=manage_acnums(acnum,Acnums2RProts)
                if proteins!=[]:
                    for prot in proteins:
                        secondary=prot
                        DatabaseGraph.is_same.create(primary,secondary)
                        j+=1

        
def clean(ObjectType):
    i=0
    for elt in ObjectType.get_all():
        ID=str(elt).split('/')[-1][:-1]
        i+=1
        logging.debug('del %s', i)
        ObjectType.delete(ID)
def getGOs(ObjectType):
    for elt in ObjectType.get_all():
        ID=str(elt).split('/')[-1][:-1]
        primary=ObjectType.get(ID)
        GODict[primary.ID]=primary
    logging.debug('GO Loaded')
        

clean(DatabaseGraph.UNIPORT)                                                                                                                                                                                                                                                                                                                                                                                                                    
getGOs(DatabaseGraph.GOTerm)
#import_GOs()
import_UNIPROTS()
