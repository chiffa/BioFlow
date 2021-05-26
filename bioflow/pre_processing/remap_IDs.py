"""
Remaps the IDs of the gene identifiers to a different organism before processing them

These are mainly interesting to for the applications where there is little information about the
genetic networks specific to the organism in question and we would like to use it as a model for
another organism (eg mice vs human) and we want to project genes associated to the phenotype in
the model organism into the networks associated to the original organism.
"""
from csv import reader as csv_reader
from csv import writer as csv_writer

from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def translate_identifiers(data_source_location, data_dump_location,
                          translation_file_location, gene_to_id_file_location):
    """
    Performs a translation of gene identifiers from one organism to another one.

    :param data_source_location: where the gene id to translate are
    :param data_dump_location: where the translated gene ids will be stored
    :param translation_file_location: where the file that contains gene id mappings is located.
        expects a .tsv file with [Source org gene id, dest org gene id, dest org gene symbol,
        confidence (1=high, 0=low), ...] per line
    :param gene_to_id_file_location: (optional) where the file that maps gene ids to HGCN names
         expects a .tsv file with [Gene stable ID, Transcript stable ID, Gene name, HGNC symbol]
    :return:
    """
    high_conf_translation_dict = {}
    low_conf_translation_dict = {}
    genes_to_ids_dict = {}

    with open(translation_file_location, 'rt') as source:
        reader = csv_reader(source, delimiter='\t')
        log.debug('org translation file header: %s' % str(next(reader)))
        for line in reader:
            if line[0] and line[1]:
                if int(line[3]):
                    # We still need to account for the confidence in mapping
                    high_conf_translation_dict[line[0]] = [line[1], line[2]]
                    # print line[0:4]
                else:
                    low_conf_translation_dict[line[0]] = [line[1], line[2]]

    high_conf_trans = []
    low_conf_trans = []

    if gene_to_id_file_location:
        with open(gene_to_id_file_location, 'rt') as source:
            reader = csv_reader(source, delimiter='\t')
            log.debug('gene to id file header: %s' % str(next(reader)))
            for line in reader:
                genes_to_ids_dict[line[2]] = line[0]

    with open(data_source_location, 'rt') as source:
        reader = csv_reader(source)
        for i, line in enumerate(reader):
            word = line[0]
            if gene_to_id_file_location:
                word = genes_to_ids_dict.get(word, 'None found')
            if word in list(high_conf_translation_dict.keys()):
                high_conf_trans.append(high_conf_translation_dict[word])
            if word in list(low_conf_translation_dict.keys()):
                low_conf_trans.append(low_conf_translation_dict[word])

    log.info("out of %s, %s were translated with high confidence,"
             " %s with low and %s were not found" % \
             (i, len(high_conf_trans), len(low_conf_trans),
              i - len(high_conf_trans) - len(low_conf_trans)))

    with open(data_dump_location, 'wt') as destination:
        writer = csv_writer(destination)
        writer.writerows((word for word in high_conf_trans))


if __name__ == "__main__":
    translation_file_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/Mouse_2_human.tsv'
    gene_to_id_file_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/mouse_look_up_table.tsv'
    data_source_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/all_significant.csv'
    data_dump_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/all_sig_hum.csv'
    translate_identifiers(data_source_location, data_dump_location,
                          translation_file_location, gene_to_id_file_location)