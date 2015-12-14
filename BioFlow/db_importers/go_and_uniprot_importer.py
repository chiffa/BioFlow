"""
This module is responsible for the import of the data parsed from the UNIPROT text file
"""
from BioFlow import main_configs as conf
from BioFlow.utils.log_behavior import logger as log
from BioFlow.bio_db_parsers.geneOntologyParser import GOTermsParser
from BioFlow.bio_db_parsers.uniprotParser import UniProtParser
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph

# Stores relations between GO IDs and the objects in the neo4j database
GO_term_memoization_dict = {}
# Stores relations between the SWISSPROT UNIPROT IDs and the neo4j database objects
Uniprot_memoization_dict = {}


def import_gene_ontology(go_terms, go_terms_structure):
    """
    Imports GOs by loading GO_Terms and GO_Terms structure from utils.GO_Structure_Parser

    :param go_terms:
    :param go_terms_structure:
    """
    # generate terms:
    go_terms, go_terms_structure = GOTermsParser().parse_go_terms(conf.GeneOntology)
    # Create Terms
    go_terms_number = len(go_terms.keys())
    for i, GO_Term in enumerate(go_terms.keys()):
        log.debug('creating GO terms: %s %%',
                  str("{0:.2f}".format(float(i) / float(go_terms_number) * 100)))
        GO_term_memoization_dict[GO_Term] = DatabaseGraph.GOTerm.create(
            ID=go_terms[GO_Term]['id'],
            Name=go_terms[GO_Term]['name'],
            displayName=go_terms[GO_Term]['name'],
            Namespace=go_terms[GO_Term]['namespace'],
            Definition=go_terms[GO_Term]['def'])

    # Create the structure between them:
    go_links_number = len(go_terms_structure)
    for i, relation in enumerate(go_terms_structure):
        log.debug('creating GO terms relations %s %%',
                  str("{0:.2f}".format(float(i) / float(go_links_number) * 100)))
        go_term_1 = GO_term_memoization_dict[relation[0]]
        go_term_2 = GO_term_memoization_dict[relation[2]]
        go_relation_type = relation[1]
        if go_relation_type == 'is_a':
            DatabaseGraph.is_a_go.create(go_term_1, go_term_2)
        if go_relation_type == 'part_of':
            DatabaseGraph.is_part_of_go.create(go_term_1, go_term_2)
        if 'regul' in go_relation_type:
            if go_relation_type == 'positively_regulates':
                DatabaseGraph.is_regulant.create(
                    go_term_1, go_term_2, controlType='ACTIVATES', ID=str(
                        'GO' + go_term_1.ID + go_term_2.ID))
            if go_relation_type == 'negatively_regulates':
                DatabaseGraph.is_regulant.create(
                    go_term_1, go_term_2, controlType='INHIBITS', ID=str(
                        'GO' + go_term_1.ID + go_term_2.ID))
            else:
                DatabaseGraph.is_regulant.create(
                    go_term_1, go_term_2, ID=str(
                        'GO' + go_term_1.ID + go_term_2.ID))


def getExistingAcnums():
    """
    Attempts to retrieve acnums from the neo4j database

    :return: dict that maps acnums to nodes in the database to which they point (Reactome proteins)

    :raise Exception: if the generator that performs a lookup of existing uniprot acnum annotation nodes is null (Reactome
    wasn't imported yet)
    """
    annot_node_with_acc_nums_generator = DatabaseGraph.AnnotNode.index.lookup(
        ptype='UniProt')
    if annot_node_with_acc_nums_generator is None:
        raise Exception(
            "Reactome was not loaded or contains no acc_num cross-references to Uniprot")
    AcnumList = {}  # acnum to AnnotNode
    for annot_node in annot_node_with_acc_nums_generator:
        if annot_node is not None:
            annot_node_ID = str(annot_node).split('/')[-1][:-1]
            AnnotObj = DatabaseGraph.vertices.get(annot_node_ID)
            AcnumList[str(AnnotObj.payload)] = AnnotObj
    if len(AcnumList) < 10:
        raise Exception(
            "Reactome was not loaded or contains no acc_num cross-references to Uniprot")
    ReactProtList = {}
    for acc_num in AcnumList.keys():
        ReactProtGen = AcnumList[acc_num].bothV()
        ReactProtList[acc_num] = []
        if ReactProtGen is not None:
            for vertex in ReactProtGen:
                if vertex is not None:
                    ReactProtList[acc_num].append(vertex)
    log.debug('1 %s', len(ReactProtList))
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
    prot_node = Uniprot_memoization_dict[CH_PROT_ID]
    annot_node = DatabaseGraph.AnnotNode.create(ptype=p_type, payload=p_load)
    DatabaseGraph.is_annotated.create(prot_node, annot_node)


def import_UNIPROTS():
    """
    Imports the whole parsed uniprot dictionary from the utils.uniprot parser into the database
    """
    Uniprot = UniProtParser(conf.up_tax_ids).parse_uniprot(conf.UNIPROT_source)
    Acnums2RProts = getExistingAcnums()
    i = 0
    j = 0
    Acnum_key_no = len(Acnums2RProts.keys()) / 100.
    UP_key_no = len(Uniprot.keys()) / 100.
    for k, CH_PROT_ID in enumerate(Uniprot.keys()):
        set1 = set(Uniprot[CH_PROT_ID]['Acnum'])
        set2 = set(Acnums2RProts.keys())
        # Create uniprot terms
        primary = DatabaseGraph.UNIPORT.create(
            ID=CH_PROT_ID,
            displayName=Uniprot[CH_PROT_ID]['Names']['Full'],
            main_connex=False)
        print CH_PROT_ID, primary
        # check inclusion in Reactome
        # TODO: if all explodes on the next import, check the line below and
        # revert the import behavior of Uniprot
        if not set1.isdisjoint(set2):
            log.debug(
                'Uniprot %s intersects Reactome on the following acnums: %s' %
                (str(CH_PROT_ID), str(set1)))
            i += len(set1)
            primary.involved = True
        # TODO: if all explodes on the next import, check the line above and
        # revert the import behavior of Uniprot
        advance_1 = i / Acnum_key_no
        advance_2 = k / UP_key_no
        log.debug(
            'loading UNIPROT: Acnums cross-linked - %.2f %% ; Total loaded: - %.2f %%' %
            (advance_1, advance_2))
        # Add the newly created uniprot to the buffer
        Uniprot_memoization_dict[CH_PROT_ID] = primary
        # Insert references to GOs
        for GO_Term in Uniprot[CH_PROT_ID]['GO']:
            if GO_Term in GO_term_memoization_dict.keys():
                secondary = GO_term_memoization_dict[GO_Term]
                DatabaseGraph.is_go_annotation.create(primary, secondary)
        # Find intersection between logged acnums and acnums pointed out by the
        # Reactome.org. Create a direct bridge between a reactome ID and a
        # SWISSPROT ID
        for acnum in Uniprot[CH_PROT_ID]['Acnum']:
            proteins = manage_acnums(acnum, Acnums2RProts)
            if proteins is not []:
                for prot in proteins:
                    secondary = prot
                    DatabaseGraph.is_same.create(primary, secondary)
                    j += 1
        for acc_num in Uniprot[CH_PROT_ID][
                'Acnum']:  # Already linked via reactome, but with a different identifier
            link_annotation(CH_PROT_ID, 'UNIPROT_Accnum', acc_num)
        link_annotation(
            CH_PROT_ID,
            'UNIPROT_Name',
            Uniprot[CH_PROT_ID]['Names']['Full'])
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
            link_annotation(
                CH_PROT_ID,
                'UNIPROT_EMBL_AC|' +
                dico['status'] +
                '|' +
                dico['type'],
                dico['Accession'].upper())
            link_annotation(
                CH_PROT_ID,
                'UNIPROT_EMBL_ID|' +
                dico['status'] +
                '|' +
                dico['type'],
                dico['ID'].upper())

        for name in Uniprot[CH_PROT_ID]['PDB']:
            link_annotation(CH_PROT_ID, 'UNIPROT_PDB', name.upper())


def getGOs():
    """
    re-loads GO
    """
    print 'recovering GO terms'
    log.debug('starting GO load')
    ObjectType = DatabaseGraph.GOTerm
    for elt in ObjectType.get_all():
        ID = str(elt).split('/')[-1][:-1]
        primary = ObjectType.get(ID)
        GO_term_memoization_dict[primary.ID] = primary
    log.debug('GO Loaded')


def getUniprots():
    """
    Pre-loads uniprots

    """
    print 'recovering UNIPROTs'
    log.debug('starting UNIPROT load')
    ObjectType = DatabaseGraph.UNIPORT
    for elt in ObjectType.get_all():
        CH_PROT_ID = elt.ID
        ID = str(elt).split('/')[-1][:-1]
        primary = ObjectType.get(ID)
        Uniprot_memoization_dict[CH_PROT_ID] = primary
    log.debug('UNIPROT Loaded')


if __name__ == "__main__":
    import_UNIPROTS()
