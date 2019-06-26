from bioflow.neo4j_db.db_io_routines import look_up_annotation_set
from csv import reader as csv_reader
from csv import writer as csv_writer


def retrieve_id_translation_table(id_list, id_ptable=None):
    if id_ptable is not None:
        translated_list = look_up_annotation_set(id_list, id_ptable)
    else:
        translated_list = look_up_annotation_set(id_list)

    output_collector = []

    for key, value in translated_list:
        output_collector.append([key, value[2], value[3]])

    return output_collector


def retrieve_base_table(path, header=False):
    output_collector = []

    with open(path, 'r') as source:
        reader = csv_reader(source)
        if header:
            reader.next()
        for line in reader:
            output_collector.append(line)

    return output_collector


def write_mapping_table(path, data_table, header=False):

    with open(path, 'w') as sink:
        writer = csv_writer(sink)
        if header:
            writer.writeline(['Source ID', 'DB ID', 'UNIPROT_ID'])
        for line in data_table:
            writer.writeline(line)


if __name__ == '__main__':
    pass