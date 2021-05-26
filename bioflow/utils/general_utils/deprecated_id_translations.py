"""
A set of large wrappers for the translation of the gene list ids to internal db ides
has been deprecated
"""
from bioflow.neo4j_db.db_io_routines import look_up_annotation_set
from csv import reader as csv_reader
from csv import writer as csv_writer
import os

# LEGACY: this is the place where the id translation from mouse to humans can go

def retrieve_id_translation_table(id_list, id_ptable=None):
    """
    Wrapper for the bioflow.neo4j_db.db_io_routines.look_up_annotation_set

    Looks up an set of annotations in the database and finds the Ids of nodes containing SWISSPROT
    proteins linked to by annotations

    :param id_list: list of payloads
    :param id_ptable:  expected type of payloads
    :return: for ids without a possible translation, instead adds an ~, ~ translation
    """
    if id_ptable is not None:
        translated_list = look_up_annotation_set(id_list, id_ptable)
    else:
        translated_list = look_up_annotation_set(id_list)

    output_collector = []

    for key, value in translated_list[1]:
        if len(value) != 0:
            output_collector.append([key, value[0][2], value[0][3]])
        else:
            output_collector.append([key, '~', '~'])

    return output_collector


def retrieve_base_table(path, header=False):
    """
    reads the table from which to translate

    :param path: where to look for the table
    :param header: if there is a header
    :return: table contents (list of lists)
    """
    output_collector = []

    with open(path, 'rt') as source:
        reader = csv_reader(source)
        if header:
            next(reader)
        for line in reader:
            output_collector.append(line)
    output_collector = [val for sublist in output_collector for val in sublist]
    return output_collector


def write_mapping_table(path, data_table, header=False):
    """
    Writes to the table where the translations are found

    :param path: where the table to write into is located
    :param data_table: contents of the data table
    :param header: contents of the header
    :return: None
    """
    with open(path, 'wt') as sink:
        writer = csv_writer(sink)
        if header:
            writer.writerow(['Source ID', 'DB ID', 'UNIPROT_ID'])
        for line in data_table:
            writer.writerow(line)


if __name__ == '__main__':
    base_table_1 = retrieve_base_table('/home/kucharav/table_1.csv')
    translated_table = retrieve_id_translation_table(base_table_1)
    write_mapping_table('/home/kucharav/translation_1.csv', translated_table)

    base_table_2 = retrieve_base_table('/home/kucharav/table_2.csv')
    write_mapping_table('/home/kucharav/translation_2.csv', retrieve_id_translation_table(base_table_2))