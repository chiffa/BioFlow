from csv import reader as csv_reader
import numpy as np
from matplotlib import pyplot as plt
import os
from itertools import combinations_with_replacement
import networkx as nx


prefix = "/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/TF 2 targets"
loc_1 = os.path.join(prefix, 'human_cellnet_grn_Apr_05_2017.csv')
loc_2 = os.path.join(prefix, 'TRRUST/trrust_rawdata.human.tsv')
marbach_prefix = os.path.join(prefix, 'Marbach2016/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks')


def open_cellnet_grn(cellnet_file):
    compendium = []
    compendium_support = []
    with open(cellnet_file, 'rb') as source:
        reader = csv_reader(source, delimiter=',')
        header = reader.next()
        # print header
        for line in reader:
            interaction_no = int(line[0])
            interaction_from = line[1]
            interaction_to = line[2]
            interaction_z_score = float(line[3])
            interaction_correlation = float(line[4])
            compendium.append((interaction_from, interaction_to))
            compendium_support.append((interaction_z_score, interaction_correlation))

    return compendium, compendium_support


def open_TRRUST(trrust_file):
    compendium = []
    compendium_support = []
    with open(trrust_file, 'rb') as source:
        reader = csv_reader(source, delimiter='\t')
        for line in reader:
            interaction_from = line[0]
            interaction_to = line[1]
            interaction_type = line[2]
            evidence = line[3].split(';')
            evidence_redundancy = len(evidence)
            compendium.append((interaction_from, interaction_to))
            compendium_support.append((evidence_redundancy))

    return compendium, compendium_support


def compare_compendiums(compendium_1, compendium_1_support, compendium_2, compendium_2_support):

    intersection = set(compendium_1).intersection(set(compendium_2))

    # print len(compendium_1), len(compendium_2)
    # print 'intersection: ', len(intersection)

    compendium_1 = np.array(compendium_1)
    compendium_2 = np.array(compendium_2)

    s1 = set(compendium_1.flatten().tolist())
    s2 = set(compendium_2.flatten().tolist())
    intersection2 = s1.intersection(s2)

    # print len(s1), len(s2)
    # print 'nodes intersection: ', len(intersection2)

    compendium_1_support = np.array(compendium_1_support).astype(np.double)
    compendium_2_support = np.array(compendium_2_support).astype(np.int)

    # plt.hist(np.log(compendium_2_support), 100)
    # plt.show()
    #
    # plt.hist(np.log(compendium_1_support[:, 0]), 100)
    # plt.show()
    #
    # plt.hist(compendium_1_support[:, 1], 100)
    # plt.show()

    return (len(s1), len(s2), len(intersection)), (len(compendium_1), len(compendium_2), len(intersection))


def short_compare_compendiums(compendium_1, compendium_2,):

    intersection = set(compendium_1).intersection(set(compendium_2))

    compendium_1 = np.array(compendium_1)
    compendium_2 = np.array(compendium_2)

    s1 = set(compendium_1.flatten().tolist())
    s2 = set(compendium_2.flatten().tolist())
    intersection2 = s1.intersection(s2)

    tf1 = set(compendium_1[:, 0].tolist())
    tf2 = set(compendium_2[:, 0].tolist())
    intersection3 = tf1.intersection(tf2)

    return len(intersection2), len(intersection), len(intersection3)


def open_marbach(marbach_file):
    compendium = []
    compendium_support = []
    with open(marbach_file, 'rb') as source:
        reader = csv_reader(source, delimiter='\t')
        for line in reader:
            interaction_from = line[0]
            interaction_to = line[1]
            if len(line) > 2:
                weight = float(line[2])
            else:
                weight = np.nan
            compendium.append((interaction_from, interaction_to))
            compendium_support.append((weight))

    return compendium, compendium_support


if __name__ == "__main__":
    c1, c1_s = open_cellnet_grn(loc_1)
    c2, c2_s = open_TRRUST(loc_2)
    # compare_compendiums(c1, c1_s, c2, c2_s)
    name_2_ci_csi = {"cellnet_grn": (c1, c1_s), "TRRUST": (c2, c2_s)}
    name_2_idx = {"cellnet_grn": 0, "TRRUST": 1}
    idx_2_name = ["cellnet_grn", "TRRUST"]
    i = 2

    for file_name in os.listdir(marbach_prefix):
        print file_name
        if file_name[0] in ['0', '1', '3']:
            shortname = file_name.split('.')[0]
            location = os.path.join(marbach_prefix, file_name)
            name_2_ci_csi[shortname] = open_marbach(location)
            name_2_idx[shortname] = i
            idx_2_name.append(shortname)
            i += 1

    size = len(idx_2_name)
    node_correlations = np.zeros((size, size))
    links_correlations = np.zeros((size, size))
    tf_correlations = np.zeros((size, size))
    for i, j in combinations_with_replacement(range(0, size), 2):
        c1 = name_2_ci_csi[idx_2_name[i]][0]
        c2 = name_2_ci_csi[idx_2_name[j]][0]
        v1, v2, v3 = short_compare_compendiums(c1, c2)
        node_correlations[i, j] = v1
        links_correlations[i, j] = v2
        tf_correlations[i, j] = v3


    print idx_2_name
    print node_correlations
    print links_correlations
    print tf_correlations


    # c, c_s = open_marbach(os.path.join(marbach_prefix, "01_neurons_fetal_brain.txt"))
    # compare_compendiums(c, c_s, c1, c1_s)