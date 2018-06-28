from csv import reader as csv_reader


origin = "/home/andrei/sources/Phosphosite/Kinase_Substrate_Dataset.tsv"


def parse_phosphosite(phoshosite_file):
    base = []
    ret_dict = {}

    with open(phoshosite_file, 'rb') as source:
        reader = csv_reader(source, delimiter='\t')
        reader.next()
        reader.next()
        reader.next()
        header = reader.next()
        counter = 0
        for line in reader:
            if line[3] == 'human' and line[8] == 'human':
                interaction_from = line[0]
                interaction_to = line[7]
                in_vivo = bool(line[-3])
                in_vitro = bool(line[-2])
                base.append(interaction_from)
                base.append(interaction_to)
                ret_dict[(interaction_from, interaction_to)] = (in_vivo, in_vitro)

    base = list(set(base))

    return ret_dict, base


if __name__ == "__main__":
    print parse_phosphosite(origin)
