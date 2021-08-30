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
                          translation_file_location, gene_to_id_file_location=None,
                          low_confidence_translations_accepted=False):
    """
    Performs a translation of gene identifiers from one organism to another one.

    The relevant data in the translation file is source stable gene id (col. 1), target stable
    gene id (col. 2), target gene name (col. 3) and confidence in translation (col. 4),
    where 0 and 1 denote respectively a low and high confidences.

    The gene to id translation is optional and will translate gene names to gene stable ids if
    the original file does not use the stable gene ids (gene names change over times and there
    are variance in their nomenclature).

    (nb: col 1 = col 0 in python)

    :param data_source_location: where the gene id to translate are
    :param data_dump_location: where the translated gene ids will be stored
    :param translation_file_location: where the file that contains gene id mappings is located.
        expects a .tsv file with [Source org gene id, dest org gene id, dest org gene symbol,
        confidence (1=high, 0=low), ...] per line
    :param gene_to_id_file_location: (optional) where the file that maps gene ids to HGCN names
         expects a .tsv file with [Gene stable ID, Transcript stable ID, Gene name, HGNC symbol].
         Adding this file will automatically trigger a translation
    :param low_confidence_translations_accepted: if set to true, the contects of lowe confidence
        mapping will be returned as well
    :return:
    """
    # The look up table format :
    # Gene stable id org 1  |   Transcript stable id org 2    |   Gene name   | HGCN symbol

    # the translation file format:
    # Gene stable ID | Human gene stable ID | Human gene name |
    # Human orthology confidence [0 low, 1 high] | Gene name | Gene description


    high_conf_translation_dict = {}  # Gene mappings with high orthology confidence
    low_conf_translation_dict = {}  # Gene mappings with low orthology confidence
    genes_to_ids_dict = {}  # Means ?

    with open(translation_file_location, 'rt') as source:
        reader = csv_reader(source, delimiter='\t')
        log.info('Parsing organism translation table')
        log.debug('org translation file: %s' % translation_file_location)
        log.debug('org translation file header: %s' % str(next(reader)))

        for line in reader:
            from_gen_id, to_gen_id, to_gene_name, confidence, from_gene_id, desc = line

            if from_gen_id and to_gen_id:  # there is a possible translation

                if int(confidence):  # 0: low, 1: high
                    high_conf_translation_dict[from_gen_id] = [to_gen_id, to_gene_name]

                else:
                    low_conf_translation_dict[from_gen_id] = [to_gen_id, to_gene_name]
        log.info('Parsing done')

    high_conf_trans = []
    low_conf_trans = []

    if gene_to_id_file_location:  # reverse translates the name of genes into stable gene ids
        with open(gene_to_id_file_location, 'rt') as source:
            reader = csv_reader(source, delimiter='\t')
            log.info('Parsing gene to id file')
            log.debug('gene to id file: %s' % translation_file_location)
            log.debug('gene to id file header: %s' % str(next(reader)))

            for line in reader:
                genes_to_ids_dict[line[2]] = line[0]

            log.info('Parsing done')

    total_lines = 0
    with open(data_source_location, 'rt') as source:
        reader = csv_reader(source, delimiter='\t')

        for i, line in enumerate(reader):
            gene_id = line[0]
            gene_weight = None
            if len(line) > 1:
                gene_weight = line[1]

            # print("attempting to translate '%s, %s'" % (gene_id, gene_weight))

            if gene_to_id_file_location:
                gene_id = genes_to_ids_dict.get(gene_id, 'None found')
                # print('gene to id translation: %s' % gene_id)

            if gene_id in list(high_conf_translation_dict.keys()):
                # print("high-confidence translation '%s'" % high_conf_translation_dict[gene_id])

                if gene_weight is not None:
                    high_conf_trans.append([high_conf_translation_dict[gene_id][1], gene_weight])
                    # print("weighted high-confidence translation '%s, %s'" %
                    #       (high_conf_translation_dict[gene_id][1], gene_weight))
                else:
                    high_conf_trans.append([high_conf_translation_dict[gene_id][1], ])
                    # print("unweighted high-confidence translation '%s'" %
                    #       high_conf_translation_dict[gene_id][1])

            if gene_id in list(low_conf_translation_dict.keys()):
                # print("high-confidence translation '%s'" % high_conf_translation_dict[gene_id])

                if gene_weight is not None:
                    low_conf_trans.append([low_conf_translation_dict[gene_id][1], gene_weight])
                    # print("weighted low-confidence translation '%s, %s'" %
                    #       (low_conf_translation_dict[gene_id][1], gene_weight))
                else:
                    low_conf_trans.append([low_conf_translation_dict[gene_id][1], ])
                    # print("unweighted low-confidence translation '%s'" %
                    #       low_conf_translation_dict[gene_id][1])

            # input('confirm line translation')

            total_lines = i

    # TODO: there are some ids that are translated both with high and low confidence, or map to
    #  the same values. Which is problematic

    log.info("out of %s ids, %s were translated with high confidence,"
             " %s with low and %s were not found" % \
             (total_lines, len(high_conf_trans), len(low_conf_trans),
              total_lines - len(high_conf_trans) - len(low_conf_trans)))

    if len(high_conf_trans) + len(low_conf_trans) < 0.2 * total_lines:
        raise Exception('Problem with translation - too few genes were translated between '
                        'organisms. Please check the format compatibility')

    with open(data_dump_location, 'wt') as destination:
        writer = csv_writer(destination, delimiter='\t')
        writer.writerows((word for word in high_conf_trans))
        if low_confidence_translations_accepted:
            writer.writerows((word for word in low_conf_trans))


if __name__ == "__main__":
    pass
