"""
Protein relationships parser
"""
# from csv import reader
from csv import reader as csv_reader
from collections import defaultdict


def parse_bio_grid(bio_grid):
    """
    Parses the given file as a BioGrid file and returns as

    :param bio_grid: the location of the BioGRID src file that needs to bprased
    :return:
    """
    ret_dict = {}
    base = []

    with open(bio_grid, 'rb') as source_file:
        biogrid_reader = csv_reader(source_file, 'excel-tab')
        biogrid_reader.next()
        for fields in biogrid_reader:
            ret_dict[tuple(fields[7:9])] = [fields[17]]
            if fields[18] != '-':
                ret_dict[tuple(fields[7:9])].append(fields[18])
            base.append(fields[7])
            base.append(fields[8])

    return ret_dict, base


def parse_hint(_hint_csv):
    """
    Reads protein-protein relationships from a HiNT database file

    :param _hint_csv: location of the HiNT database tsv file
    :return: {UP_Identifier:[UP_ID1, UP_ID2, ...]}
    """
    local_relations = defaultdict(list)

    with open(_hint_csv, 'r') as source_file:
        hint_reader = csv_reader(source_file, delimiter='\t')
        hint_reader.next()
        for i, fields in enumerate(hint_reader):
            if fields[2] != fields[3]:
                local_relations[fields[3]].append(fields[2])
                local_relations[fields[2]].append(fields[3])
    return dict(local_relations)
