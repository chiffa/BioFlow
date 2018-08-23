"""
This module is responsible for the import of the data parsed from the UNIPROT text file
"""
from bioflow import main_configs
from bioflow.utils.log_behavior import get_logger
from bioflow.bio_db_parsers.geneOntologyParser import GOTermsParser
from bioflow.bio_db_parsers.uniprotParser import UniProtParser
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.neo4j_db.db_io_routines import get_db_id
from bioflow.internal_configs import neo4j_names_dict

log = get_logger(__name__)

# TODO: REFACTORING: this should be refactored into a class to wrap memoization dicts
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
    # Create Terms
    go_terms_number = len(go_terms.keys())
    log.info('Starting to importing %s GO terms' % go_terms_number)

    for i, (GO_Term, term) in enumerate(go_terms.iteritems()):
        log.debug('creating GO terms: %s %%',
                  str("{0:.2f}".format(float(i) / float(go_terms_number) * 100)))
        if i*20 % go_terms_number < 20:
            log.info('GO terms import: %s %%',
                     "{0:.2f}".format(float(i) / float(go_terms_number) * 100))
        GO_term_memoization_dict[GO_Term] = DatabaseGraph.create(neo4j_names_dict['GO Term'],
                             {'legacyId': term['id'],
                              'Name': term['name'],
                              'displayName': term['name'],
                              'Namespace': term['namespace'],
                              'Definition': term['def']})

    # Create the structure between them:
    go_links_number = len(go_terms_structure)
    log.info('Starting to import %s GO terms links' % go_links_number)

    for i, relation in enumerate(go_terms_structure):
        log.debug('creating GO terms relations %s %%',
                  str("{0:.2f}".format(float(i) / float(go_links_number) * 100)))
        if i*20 % go_links_number < 20:
            log.info('GO terms linking: %s %%',
                     "{0:.2f}".format(float(i) / float(go_links_number) * 100))
        go_term_1 = GO_term_memoization_dict[relation[0]]
        go_term_2 = GO_term_memoization_dict[relation[2]]
        go_relation_type = relation[1]

        if go_relation_type == 'is_a':
            DatabaseGraph.link(get_db_id(go_term_1), get_db_id(go_term_2), 'is_a_go')

        if go_relation_type == 'part_of':
            DatabaseGraph.link(get_db_id(go_term_1), get_db_id(go_term_2), 'is_part_of_go')

        if 'regul' in go_relation_type:
            if go_relation_type == 'positively_regulates':
                DatabaseGraph.link(get_db_id(go_term_1), get_db_id(go_term_2),
                                   'is_regulant',
                                   {'controlType': 'ACTIVATES',
                                    'ID': str('GO' + go_term_1.properties['legacyId'] +
                                              go_term_2.properties['legacyId'])})

            if go_relation_type == 'negatively_regulates':
                DatabaseGraph.link(get_db_id(go_term_1), get_db_id(go_term_2),
                                   'is_regulant',
                                   {'controlType': 'INHIBITS',
                                    'ID': str('GO' + go_term_1.properties['legacyId'] +
                                              go_term_2.properties['legacyId'])})

            else:
                DatabaseGraph.link(get_db_id(go_term_1), get_db_id(go_term_2),
                                   'is_regulant',
                                   {'ID': str('GO' + go_term_1.properties['legacyId'] +
                                              go_term_2.properties['legacyId'])})


def pull_up_acc_nums_from_reactome():
    """
    Attempts to retrieve accession numbers nums from the neo4j database

    :return: dict that maps acnums to nodes in the database to which they point (Reactome proteins)
    :raise Exception: if the generator that performs a lookup of existing uniprot acnum annotation
    nodes is null (Reactome wasn't imported yet) or a wrong uniprot is crosslinked with a wrong
    reactome
    """
    acc_num_annot_nodes = DatabaseGraph.find({'type': 'UniProt'}, 'Annotation')

    if acc_num_annot_nodes is None:
        raise Exception(
            "Reactome was not loaded or contains no acc_num cross-references to Uniprot")

    reactome_proteins = {}

    breakpoints = 300
    total_nodes = len(acc_num_annot_nodes)

    for i, annotation_node in enumerate(acc_num_annot_nodes):

        if i % breakpoints:
            log.info("\t cross-linking %.2f %% complete" % (float(i)/float(total_nodes)*100.))

        tag = annotation_node.properties['tag']
        node = DatabaseGraph.get_from_annotation_tag(tag, 'UniProt')
        reactome_proteins[tag] = node

        if '-' in tag:
            reactome_proteins[tag.split('-')[0]] = node
        if '_' in tag:
            reactome_proteins[tag.split('_')[0]] = node
        if '.' in tag:
            reactome_proteins[tag.split('.')[0]] = node
        if ',' in tag:
            reactome_proteins[tag.split(',')[0]] = node

        if reactome_proteins[tag] != 1:
            log.debug('Cross-linking reactome v.s. acc_num %s mapped to %s proteins',
                      tag, len(reactome_proteins[tag]))

    log.info('Cross-linked %s proteins from reactome v.s. Uniprot', len(reactome_proteins))

    return reactome_proteins


def manage_acc_nums(acc_num, acc_num_2_reactome_proteins):
    """
    implements an accelerated dict-assisted buffering on retrieval of relations between accession
     numbers and uniprot IDs

    :param acc_num: accession number
    :param acc_num_2_reactome_proteins: Dictionary mapping accession numbers to protein IDs
    :return: protein list if acession number has been buffered, nothing otherwise
    """
    if acc_num not in acc_num_2_reactome_proteins.keys():
        return []
    if acc_num in acc_num_2_reactome_proteins.keys():
        return acc_num_2_reactome_proteins[acc_num]


def link_annotation(uniprot_id, p_type, p_load, preferential=False):
    """
    Links a uniprot node to an annotation node

    :param uniprot_id: uniprot/swissprot ID
    :param p_type: type of the annotation
    :param p_load: content of the annotation
    """
    prot_node = Uniprot_memoization_dict[uniprot_id]
    DatabaseGraph.attach_annotation_tag(prot_node.id, p_load, p_type, preferential)


def insert_uniprot_annotations(swiss_prot_id, data_container):
    """
    Inserts all the annotations for a single swiss_prot_id contained in a supplied container

    :param swiss_prot_id:
    :param data_container:
    :return:
    """
    # TODO: single commit annotation insertion would be faster => refactor
    link_annotation(swiss_prot_id, 'UNIPROT_Name',
                    data_container['Names']['Full'], preferential=True)

    for acc_num in data_container['Acnum']:
        link_annotation(swiss_prot_id, 'UNIPROT_Accnum', acc_num)

    for name in data_container['Names']['AltNames']:
        link_annotation(swiss_prot_id, 'UNIPROT_AltName', name.upper())

    for name in data_container['GeneRefs']['Names']:
        link_annotation(swiss_prot_id, 'UNIPROT_GeneName', name.upper(), preferential=True)

    for name in data_container['GeneRefs']['AltNames']:
        link_annotation(swiss_prot_id, 'UNIPROT_AltGeneName', name.upper())

    for name in data_container['GeneRefs']['OrderedLocusNames']:
        link_annotation(swiss_prot_id, 'UNIPROT_GeneOL', name.upper())

    for name in data_container['GeneRefs']['ORFNames']:
        link_annotation(swiss_prot_id, 'UNIPROT_GeneORF', name.upper())

    for name in set(data_container['Ensembl']):
        link_annotation(swiss_prot_id, 'UNIPROT_Ensembl', name.upper())

    for dico in data_container['EMBL']:
        link_annotation(
            swiss_prot_id,
            'UNIPROT_EMBL_AC|%s|%s' % (dico['status'], dico['type']),
            dico['Accession'].upper())
        link_annotation(
            swiss_prot_id,
            'UNIPROT_EMBL_ID||%s|%s' % (dico['status'], dico['type']),
            dico['ID'].upper())

    for name in data_container['PDB']:
        link_annotation(swiss_prot_id, 'UNIPROT_PDB', name.upper())


def import_uniprots(uniprot, reactome_acnum_bindings):
    """
    Imports the whole parsed uniprot dictionary from the utils.uniprot parser into the database

    :param uniprot:
    :param reactome_acnum_bindings:
    """

    cross_links = 0
    is_same_links_no = 0
    acc_nums_no = len(reactome_acnum_bindings.keys()) / 100.
    up_no = float(len(uniprot.keys()))

    log.info('Starting to import %s uniprot nodes' % up_no)

    for sp_id_num, (swiss_prot_id, data_container) in enumerate(uniprot.iteritems()):

        set1 = set(data_container['Acnum'])
        set2 = set(reactome_acnum_bindings.keys())
        # Create uniprot terms

        uniprot_node = DatabaseGraph.create(
            neo4j_names_dict['UNIPROT'],
            {'legacyId': swiss_prot_id,
             'displayName': data_container['Names']['Full'],
             'main_connex': False})

        log.debug('uniprot %s created @ %s', swiss_prot_id, uniprot_node)

        if not set1.isdisjoint(set2):
            log.debug('Uniprot %s intersects Reactome on the following acnums: %s',
                      str(swiss_prot_id), str(set1))
            cross_links += len(set1)
            uniprot_node.involved = True

        log.debug(
            'loading UNIPROT: Acnums cross-linked - %.2f %% ; Total loaded: - %.2f %%',
            cross_links / acc_nums_no, sp_id_num / up_no *100)

        if sp_id_num*20 % up_no < 20:
            log.info('Uniprots nodes load: %s %%', sp_id_num / up_no * 100 )

        # Add the newly created uniprot to the buffer
        Uniprot_memoization_dict[swiss_prot_id] = uniprot_node

        # Insert references to GOs
        for GO_Term in data_container['GO']:
            if GO_Term in GO_term_memoization_dict.keys():
                linked_go_term = GO_term_memoization_dict[GO_Term]
                DatabaseGraph.link(get_db_id(uniprot_node), get_db_id(linked_go_term), 'is_go_annotation')

        for acnum in data_container['Acnum']:
            proteins = manage_acc_nums(acnum, reactome_acnum_bindings)
            if proteins is not []:
                for prot in proteins:
                    secondary = prot
                    DatabaseGraph.link(get_db_id(uniprot_node), get_db_id(secondary), 'is_same')
                    is_same_links_no += 1

        insert_uniprot_annotations(swiss_prot_id, data_container)


def memoize_go_terms():
    """
    loads go terms from the
    """
    for node in DatabaseGraph.get_all(neo4j_names_dict['GO Term']):
        GO_term_memoization_dict[node.properties['legacyId']] = node


def memoize_uniprots():
    """
    Pre-loads uniprots
    """
    for node in DatabaseGraph.get_all(neo4j_names_dict['UNIPROT']):
        Uniprot_memoization_dict[node.properties['legacyId']] = node


if __name__ == "__main__":
    _uniprot = UniProtParser(main_configs.up_tax_ids).parse_uniprot(main_configs.uniprot_path)
    _reactome_acnum_bindings = pull_up_acc_nums_from_reactome()
    import_uniprots(_uniprot, _reactome_acnum_bindings)
