__author__ = 'ank'

from PolyPharma.configs import RNA_source
import csv

def read_RNA_seq_counts():
    retlist = []
    print 'opening: ', RNA_source
    with open(RNA_source,'rb') as source_file:
        source_reader = csv.reader(source_file, dialect='excel-tab')
        for row in source_reader:
            retlist.append(row)
    return retlist

def filter_id_only(parse_table):
    return [sublist[0] for sublist in parse_table if 'ENS' in sublist[0]]

if __name__ == "__main__":
    print filter_id_only(read_RNA_seq_counts())