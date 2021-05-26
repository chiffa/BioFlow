from csv import reader as csv_reader
import numpy as np
from collections import defaultdict
import os
import bioflow.configs.main_configs as mc


def parse_TRRUST(trrust_file):
    base = []
    ret_dict = {}

    with open(trrust_file, 'rt') as source:
        reader = csv_reader(source, delimiter='\t')
        for line in reader:
            interaction_from = line[0]
            interaction_to = line[1]
            interaction_type = line[2]
            evidence = line[3].split(';')
            evidence_redundancy = len(evidence)
            base.append(interaction_to)
            base.append(interaction_from)
            ret_dict[(interaction_from, interaction_to)] = evidence_redundancy

    base = list(set(base))

    return ret_dict, base


def parse_cellnet_grn(cellnet_file):
    base = []
    ret_dict = {}

    with open(cellnet_file, 'rt') as source:
        reader = csv_reader(source, delimiter=',')
        header = next(reader)
        # print header
        for line in reader:
            interaction_no = int(line[0])
            interaction_from = line[1]
            interaction_to = line[2]
            interaction_z_score = float(line[3])
            interaction_correlation = float(line[4])
            base.append(interaction_to)
            base.append(interaction_from)
            ret_dict[(interaction_from, interaction_to)] = interaction_correlation

    base = list(set(base))

    return ret_dict, base


def parse_marbach(marbach_prefix, parse_mode='mean'):
    """


    :param marbach_prefix:
    :param parse_mode: ['mean', 'one', 'all']
    :raise Exception: unsupported parse mode
    :return:
    """

    def open_marbach(marbach_file, insertion_index):
        with open(marbach_file, 'rt') as source:
            reader = csv_reader(source, delimiter='\t')
            for line in reader:
                interaction_from = line[0]
                interaction_to = line[1]

                if len(line) > 2:
                    weight = np.abs(float(line[2]))
                else:
                    weight = np.nan

                master_accumulator[(interaction_from, interaction_to)][insertion_index] = weight


    master_accumulator = defaultdict(lambda: np.zeros((32, )))

    for file_name in os.listdir(marbach_prefix):
        # print(file_name)
        if file_name[0] in ['0', '1', '3'] and file_name[-4:] == ".txt":
            location = os.path.join(marbach_prefix, file_name)
            open_marbach(location, insertion_index=int(file_name[:2])-1)

    ret_dict = {}
    base = []

    for pair, weight_array in master_accumulator.items():
        if parse_mode == 'mean':
            summary = sum(weight_array != 0) >= 16

        elif parse_mode == 'one':
            summary = any(weight_array != 0)

        elif parse_mode == 'all':
            summary = all(weight_array != 0)

        else:
            raise Exception("unsupported parse mode: %s" % parse_mode)

        if summary:
            ret_dict[pair] = np.mean(weight_array)
            base.append(pair[0])
            base.append(pair[1])

    base = list(set(base))

    return ret_dict, base
