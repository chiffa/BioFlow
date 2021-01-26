from csv import reader as csv_reader


def parse_phosphosite(phoshosite_file, organism):
    """
    Parses the phosphocite tsv file

    :param phoshosite_file:
    :param organism:
    :return:
    """
    base = []
    ret_dict = {}

    with open(phoshosite_file, 'rt', encoding='latin-1') as source:  # TODO: formalize the encoding
        # usage
        reader = csv_reader(source, delimiter='\t')
        next(reader)
        next(reader)
        next(reader)
        header = next(reader)
        for line in reader:
            if line[3] == organism and line[8] == organism:
                interaction_from = line[0]
                interaction_to = line[7]
                in_vivo = bool(line[-3])
                in_vitro = bool(line[-2])
                base.append(interaction_from)
                base.append(interaction_to)
                ret_dict[(interaction_from, interaction_to)] = (in_vivo, in_vitro)

    base = list(set(base))

    return ret_dict, base
