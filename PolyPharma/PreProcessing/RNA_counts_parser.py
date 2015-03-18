"""
Module responsible for importing raw RNA-seq data and running tests on it to get statistically significant differences
in genes.
"""
__author__ = 'ank'
from PolyPharma.configs2 import RNA_source, Dumps, Outputs
from PolyPharma.neo4j_analyzer.DB_IO_Routines import look_up_Annot_Node, recover_annotation
from PolyPharma.neo4j_analyzer.IO_Routines import dump_object
from collections import defaultdict
import numpy as np
from scipy import special
from csv import reader, writer
from pprint import PrettyPrinter


def import_counts_table(counts_size):
    """
    Imports the counts tables from a csv file. The csv file has to contain gene identifiers, Uxon lengths and counts
    for each identifier and each experiment. Line is gene, column is experiment

    :param counts_size: number of experiments we want to account for
    """
    gene_names = np.zeros((1, 2))
    uxon_lenghts = np.zeros((1, 1))
    table = np.zeros((counts_size, 1)).T
    initset = set()

    print RNA_source
    with open(RNA_source, 'rb') as source_file:
        rdr = reader(source_file, 'excel-tab')
        rdr.next()
        for row in rdr:
            DB_ID = look_up_Annot_Node(row[0])
            if DB_ID:
                DB_ID = DB_ID[0][2]
                gene_names = np.concatenate((gene_names, np.array([[row[0], DB_ID]])))
                initset.add(DB_ID)
                rrw = [float(f) for f in row[2 : counts_size + 2] ]
                table = np.concatenate( (table, np.array([rrw], dtype = np.float64)))
                uxon_lenghts = np.concatenate( (uxon_lenghts, np.array([[ float(row[1])]])))
    print table.shape
    return gene_names[1: , :], uxon_lenghts[1: , :], table[1: , :]


def write_out_annotation(Node_Id_list):

    with open(Outputs.cross_refs,'w') as outfile:
        _writer = writer(outfile, dialect='excel-tab')
        _writer.writerows(Node_Id_list)


def counts_filter(table, experiment_groups, filter_level):
    """
    Filters out genes for whom the count levels are too low in all the samples


    :param table: table of counts; line is gene, row is expriment
    :param experiment_groups: groups expriments into repeats
    :param filter_level: minimum amount of counts in any experience for which we are going to accept a gene
    :return: a boolean array which is True if Gene passed that test
    :rtype : np.array of dtype bool
    """
    mins = np.zeros((table.shape[0], len(experiment_groups)))

    for i, group in enumerate(experiment_groups):
        mins[:, i] = np.min(table[:, group], axis=1)
    filtr = np.max(mins, axis=1) > filter_level -1

    return filtr


def translate_to_rpkm(uxon_length, table):
    """
    Translates counts to the RPKMs

    :param uxon_length: vector of gene lengths
    :param table: table of counts
    :return: table of RPKMS
    """
    total_series_reads = np.sum(table, axis=0)
    re_table = table / total_series_reads[np.newaxis, :] / uxon_length[:, 0][:, np.newaxis] * 10e12
    return re_table


def erf_test(rpkm_table, experiment_groups, intergorups, target_p_value=0.05):
    """
    performs a test that uses the errfunction to determine if we can reject the hypothesis that all the genes
    are sampled from the same distribution

    :param rpkm_table: table of the rpkm values
    :param experiment_groups: groups on indexes
    :param intergorups: the groups between which we want to do the comparisons
    :param target_p_value: p_value with which we want to be albe to reject the null hypothesis
    """
    groups_means = np.zeros((rpkm_table.shape[0], len(experiment_groups)))
    groups_var = np.zeros((rpkm_table.shape[0], len(experiment_groups)))

    for i, group in enumerate(experiment_groups):
        groups_means[:, i] = np.mean(rpkm_table[:, group], axis=1)
        groups_var[:, i] = np.var(rpkm_table[:, group], axis=1) / estimator_dilatation_table[len(group)]**2

    group_comparison = []
    for intergroup in intergorups:
        ig_mean = np.fabs(groups_means[:, intergroup[0]] - groups_means[:, intergroup[1]])
        ig_std = np.sqrt(groups_var[:, intergroup[0]] + groups_var[:, intergroup[1]])
        p_val = 1-special.erf(ig_mean / ig_std)
        sorted_p_vals = np.sort(p_val, axis=0)
        li = np.array(range(0, sorted_p_vals.shape[0])) * target_p_value / sorted_p_vals.shape[0]
        pre_filtr = sorted_p_vals <= li
        refined_threshold =  np.max(sorted_p_vals[pre_filtr])
        filtr = p_val < refined_threshold
        group_comparison.append((p_val, filtr))

    return  group_comparison,


def run_test_suite(experiments, experimental_groups, intergroups, count_filter_level, false_discovery_rate):
    """
    Performs the full test suite, import and export included

    :param experiments: number of experiments
    :param experimental_groups: experiment groupings
    :param intergroups: groups to be compared
    :param count_filter_level: minimum counts to run statistics
    :param false_discovery_rate: desired false discovery rate
    """
    names, lengths, counts = import_counts_table(experiments)

    flter = counts_filter(counts, experimental_groups,  filter_level=count_filter_level)

    names = names[flter, :]
    lengths = lengths[flter, :]
    counts = counts[flter, :]

    rpkms = translate_to_rpkm(lengths, counts)

    testres = erf_test(rpkms, experimental_groups, intergroups, false_discovery_rate)

    filtr1, filtr2 = testres[0][1], testres[1][1]

    dump_object(Dumps.RNA_seq_counts_compare, [names[filtr1, 1], names[filtr2, 1]])

    return


if __name__ == "__main__":

    pp = PrettyPrinter(indent=4)

    pre_dict = {    1:0.80,
                    2:0.89,
                    3:0.92,
                    4:0.94,
                    5:0.95,
                    6:0.96,
                    7:0.965,
                    8:0.97
                }

    estimator_dilatation_table  = defaultdict(lambda:1)
    estimator_dilatation_table.update(pre_dict)


    exp_groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    intergroups = [[0, 1], [0, 2]]

    run_test_suite(9, exp_groups, intergroups, 10, 0.05)

    #
    # names, lengths, counts = import_counts_table(9)
    # flter = counts_filter(counts, exp_groups,  filter_level=10)
    #
    # names = names[flter, :]
    # lengths = lengths[flter, :]
    # counts = counts[flter, :]
    #
    # rpkms = translate_to_rpkm(lengths, counts)
    #
    # annot = np.array(recover_annotation(transcription,'UNIPROT_Ensembl')[1])
    #
    # fltr2 = np.in1d(names[:,0], annot[:,0])
    # fltr3 = np.in1d(annot[:,0], names[:,0])
    #
    # reslist = np.concatenate((names[fltr2, :], rpkms[fltr2, :]), axis=1)
    # reannot = annot[fltr3,:]
    #
    #
    # write_out_annotation(np.concatenate((reannot[reannot[:,0].argsort()],reslist[reslist[:,0].argsort()]), axis=1).tolist())
    #
    #

