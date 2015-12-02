"""
Set of tools to work with HiNT database
Created on Jul 10, 2013
Last updated on Nov 6, 2015

:author: Andrei Kucharavy
"""

from csv import reader as csv_reader
from collections import defaultdict
from BioFlow.main_configs import Hint_csv
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph


# TODO HiNT is subject to modification and thus parsing should be tested
# and shown
def get_prot2prot_rels(_hint_csv):
    """
    Reads protein-protein relationships from a HiNT database file

    :param _hint_csv: location of the HiNT database tsv file
    :return: {UP_Identifier:[UP_ID1, UP_ID2, ...]}
    """
    local_relations = defaultdict(list)
    with open(_hint_csv, 'r') as source_docu:
        rdr = csv_reader(source_docu, delimiter='\t')
        rdr.next()
        for i, fields in enumerate(rdr):
            if fields[2] != fields[3]:
                local_relations[fields[3]].append(fields[2])
                local_relations[fields[2]].append(fields[3])
    return dict(local_relations)


def get_uniprots():  # TODO: refactoring requires moving to the Database IO module
    """
    Connects to the Graph database and pulls out all of the uniprots by their identifiers

    :return:
    """
    uniprot_dict = {}
    for elt in DatabaseGraph.UNIPORT.get_all():
        ID = str(elt).split('/')[-1][:-1]
        print ID, elt._id  # TODO: check for potential simplification
        primary = DatabaseGraph.UNIPORT.get(ID)
        uniprot_dict[str(primary.ID).split('_')[0]] = primary
    return uniprot_dict


# TODO: dissociate indentification of elements that need to be joined and
# graph database action
def cross_ref_HiNT(flush):
    """
    Pulls Hint relationships and connects Uniprots in the database

    :param flush: if True, relationships are pushed to the actual graph database
    :return:
    """
    RelationDict = get_prot2prot_rels(Hint_csv)
    UniProtRefDict = get_uniprots()
    Treated = set()
    i = 0
    for key in UniProtRefDict.keys():
        if key in RelationDict.keys():
            Treated.add(key)
            for subkey in RelationDict[key]:
                if subkey in UniProtRefDict.keys() and subkey not in Treated:
                    i += 1
                    print 'HINT links:', key, subkey
                    if flush:
                        DatabaseGraph.is_interacting.create(
                            UniProtRefDict[key], UniProtRefDict[subkey])
    print i, len(Treated)
