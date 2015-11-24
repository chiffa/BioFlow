'''
Created on Jul 5, 2013

@author: andrei
'''
import logging

from BioFlow.data_parsers.GO_Structure_Parser import fill_GO_Terms
from BioFlow.data_parsers.UNIPROT_Parser import parse_uniprot
from BioFlow.neo4j_Declarations.Graph_Declarator import DatabaseGraph

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


GODict = {} # Stores relations between GO IDs and the objects in the neo4j database
UniprotDict = {} # Stores relations between the SWISSPROT UNIPROT IDs and the neo4j database objects


def import_GOs():
    """
    Imports GOs by loading GO_Terms and GO_Terms structure from Utils.GO_Structure_Parser
    """
    # generate terms:
    GO_Terms, GO_Terms_Structure = fill_GO_Terms()
    # Create Terms
    leng = len(GO_Terms.keys())
    i = 0
    for GO_Term in GO_Terms.keys():
        i += 1
        logging.debug('GO %s %%', str("{0:.2f}".format(float(i) / float(leng) * 100)))
        primary = DatabaseGraph.GOTerm.create(ID = GO_Terms[GO_Term]['id'],
                                              Name = GO_Terms[GO_Term]['name'],
                                              displayName = GO_Terms[GO_Term]['name'],
                                              Namespace = GO_Terms[GO_Term]['namespace'],
                                              Definition = GO_Terms[GO_Term]['def'])
        GODict[GO_Term] = primary
    # Create the structure between them:
    leng = len(GO_Terms_Structure)
    i = 0
    for relation in GO_Terms_Structure:
        i += 1
        logging.debug('inter GO relations %s %%', str("{0:.2f}".format(float(i) / float(leng) * 100)))
        primary = GODict[relation[0]]
        secondary = GODict[relation[2]]
        Type = relation[1]
        if Type == 'is_a':
            DatabaseGraph.is_a_go.create(primary, secondary)
        if Type == 'part_of':
            DatabaseGraph.is_part_of_go.create(primary, secondary)
        if 'regul' in Type:
            if Type == 'positively_regulates':
                DatabaseGraph.is_regulant.create(primary, secondary,
                                                 controlType = 'ACTIVATES',
                                                 ID = str('GO' + primary.ID + secondary.ID))
            if Type == 'negatively_regulates':
                DatabaseGraph.is_regulant.create(primary, secondary,
                                                 controlType = 'INHIBITS',
                                                 ID = str('GO' + primary.ID + secondary.ID))
            else:
                DatabaseGraph.is_regulant.create(primary, secondary,
                                                 ID = str('GO' + primary.ID + secondary.ID))


def getExistingAcnums():
    """
    Attempts to retrieve acnums from the neo4j database

    :return: dict that maps acnums to nodes in the database to which they point (Reactome proteins)

    :raise Exception: if the generator that performs a lookup of existing uniprot acnum annotation nodes is null (Reactome
    wasn't imported yet)
    """
    annot_node_with_acc_nums_generator = DatabaseGraph.AnnotNode.index.lookup(ptype = 'UniProt')
    if annot_node_with_acc_nums_generator is None:
        raise Exception("Reactome was not loaded or contains no acc_num cross-references to Uniprot")
    AcnumList = {} #acnum to AnnotNode
    for annot_node in annot_node_with_acc_nums_generator:
        if annot_node is not None:
            annot_node_ID = str(annot_node).split('/')[-1][:-1]
            AnnotObj = DatabaseGraph.vertices.get(annot_node_ID)
            AcnumList[str(AnnotObj.payload)] = AnnotObj
    if len(AcnumList) < 10:
        raise Exception("Reactome was not loaded or contains no acc_num cross-references to Uniprot")
    ReactProtList = {}
    for acc_num in AcnumList.keys():
        ReactProtGen = AcnumList[acc_num].bothV()
        ReactProtList[acc_num] = []
        if ReactProtGen is not None:
            for vertex in ReactProtGen:
                if vertex is not None:
                    ReactProtList[acc_num].append(vertex)
    logging.debug('1 %s', len(ReactProtList))
    return ReactProtList
    

def manage_acnums(acnum, Acnum2RProts):
    """
    implements an accelerated dict-assisted buffering on retrieval of relations between accession numbers and uniprot IDs

    :param acnum: accession number
    :param Acnum2RProts: Buffering dictionary mapping accession numbers to protein IDs
    :return: protein list if acession number has been buffered, nothing otherwise
    """
    if acnum not in Acnum2RProts.keys():
        return []
    if acnum in Acnum2RProts.keys():
        return Acnum2RProts[acnum]


def link_annotation(CH_PROT_ID, p_type, p_load):
    """
    Links a uniprote node to an annotation node

    :param CH_PROT_ID: swissprot ID
    :param p_type: type of the annotation
    :param p_load: content of the annotation
    """
    prot_node = UniprotDict[CH_PROT_ID]
    annot_node = DatabaseGraph.AnnotNode.create(ptype = p_type, payload = p_load)
    DatabaseGraph.is_annotated.create(prot_node, annot_node)


def import_UNIPROTS():
    """
    Imports the whole parsed uniprot dictionary from the utils.uniprot parser into the database
    """
    Uniprot = parse_uniprot()
    Acnums2RProts = getExistingAcnums()
    i = 0
    j = 0
    Acnum_key_no = len(Acnums2RProts.keys())/100.
    UP_key_no = len(Uniprot.keys())/100.
    for k, CH_PROT_ID in enumerate(Uniprot.keys()):
        set1 = set(Uniprot[CH_PROT_ID]['Acnum'])
        set2 = set(Acnums2RProts.keys())
        #Create uniprot terms
        primary = DatabaseGraph.UNIPORT.create(ID = CH_PROT_ID,
                                       displayName = Uniprot[CH_PROT_ID]['Names']['Full'],
                                       main_connex = False)
        print CH_PROT_ID, primary
        # check inclusion in Reactome
        # TODO: if all explodes on the next import, check the line below and revert the import behavior of Uniprot
        if not set1.isdisjoint(set2):
            logging.debug('Uniprot %s intersects Reactome on the following acnums: %s'%(str(CH_PROT_ID), str(set1)))
            i += len(set1)
            primary.involved = True
        # TODO: if all explodes on the next import, check the line above and revert the import behavior of Uniprot
        advance_1 = i / Acnum_key_no
        advance_2 = k / UP_key_no
        logging.debug('loading UNIPROT: Acnums cross-linked - %.2f %% ; Total loaded: - %.2f %%' % (advance_1, advance_2 ))
        # Add the newly created uniprot to the buffer
        UniprotDict[CH_PROT_ID] = primary
        # Insert references to GOs
        for GO_Term in Uniprot[CH_PROT_ID]['GO']:
            if GO_Term in GODict.keys():
                secondary = GODict[GO_Term]
                DatabaseGraph.is_go_annotation.create(primary, secondary)
        # Find intersection between logged acnums and acnums pointed out by the Reactome.org. Create a direct bridge between a reactome ID and a SWISSPROT ID
        for acnum in Uniprot[CH_PROT_ID]['Acnum']:
            proteins = manage_acnums(acnum, Acnums2RProts)
            if proteins is not []:
                for prot in proteins:
                    secondary = prot
                    DatabaseGraph.is_same.create(primary, secondary)
                    j += 1
        for acc_num in Uniprot[CH_PROT_ID]['Acnum']: # Already linked via reactome, but with a different identifier
            link_annotation(CH_PROT_ID, 'UNIPROT_Accnum', acc_num)
        link_annotation(CH_PROT_ID, 'UNIPROT_Name', Uniprot[CH_PROT_ID]['Names']['Full'])
        for name in Uniprot[CH_PROT_ID]['Names']['AltNames']:
            link_annotation(CH_PROT_ID, 'UNIPROT_Name', name.upper())

        for name in Uniprot[CH_PROT_ID]['GeneRefs']['Names']:
            link_annotation(CH_PROT_ID, 'UNIPROT_GeneName', name.upper())
        for name in Uniprot[CH_PROT_ID]['GeneRefs']['OrderedLocusNames']:
            link_annotation(CH_PROT_ID, 'UNIPROT_GeneOL', name.upper())
        for name in Uniprot[CH_PROT_ID]['GeneRefs']['ORFNames']:
            link_annotation(CH_PROT_ID, 'UNIPROT_GeneORF', name.upper())

        for name in set(Uniprot[CH_PROT_ID]['Ensembl']):
            link_annotation(CH_PROT_ID, 'UNIPROT_Ensembl', name.upper())

        for dico in Uniprot[CH_PROT_ID]['EMBL']:
            link_annotation(CH_PROT_ID, 'UNIPROT_EMBL_AC|'+dico['status']+'|'+dico['type'], dico['Accession'].upper())
            link_annotation(CH_PROT_ID, 'UNIPROT_EMBL_ID|'+dico['status']+'|'+dico['type'], dico['ID'].upper())

        for name in Uniprot[CH_PROT_ID]['PDB']:
            link_annotation(CH_PROT_ID, 'UNIPROT_PDB', name.upper())



def getGOs():
    """
    re-loads GO
    """
    print 'recovering GO terms'
    logging.debug('starting GO load')
    ObjectType = DatabaseGraph.GOTerm
    for elt in ObjectType.get_all():
        ID = str(elt).split('/')[-1][:-1]
        primary = ObjectType.get(ID)
        GODict[primary.ID] = primary
    logging.debug('GO Loaded')


def getUniprots():
    """
    Pre-loads uniprots

    """
    print 'recovering UNIPROTs'
    logging.debug('starting UNIPROT load')
    ObjectType = DatabaseGraph.UNIPORT
    for elt in ObjectType.get_all():
        CH_PROT_ID = elt.ID
        ID = str(elt).split('/')[-1][:-1]
        primary = ObjectType.get(ID)
        UniprotDict[CH_PROT_ID] = primary
    logging.debug('UNIPROT Loaded')
    

if __name__ == "__main__":
    import_UNIPROTS()