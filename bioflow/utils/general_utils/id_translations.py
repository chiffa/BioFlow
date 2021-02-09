from bioflow.neo4j_db.db_io_routines import look_up_annotation_set
from bioflow.configs.bioflow_home import internal_storage
from csv import reader as csv_reader
from csv import writer as csv_writer
import os


def retrieve_id_translation_table(id_list, id_ptable=None):
    if id_ptable is not None:
        translated_list = look_up_annotation_set(id_list, id_ptable)
    else:
        translated_list = look_up_annotation_set(id_list)

    output_collector = []

    print(translated_list[1])

    for key, value in translated_list[1]:
        # print key
        # print value
        if len(value) != 0:
            # print value[0][2], value[0][3]
            output_collector.append([key, value[0][2], value[0][3]])
        else:
            output_collector.append([key, '~', '~'])

    return output_collector


def retrieve_base_table(path, header=False):
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