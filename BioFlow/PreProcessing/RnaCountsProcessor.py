"""
Module responsible for importing raw RNA-seq per-gene counts and pulling out statistically
significantly different genes
"""
from BioFlow.configs2 import RNA_source, Dumps
from BioFlow.neo4j_analyzer import DB_IO_Routines
from BioFlow.neo4j_analyzer.IO_Routines import dump_object
from collections import defaultdict
import numpy as np
from scipy.stats import t
from csv import reader
from pprint import PrettyPrinter


pre_dict = {1: 0.80,
            2: 0.89,
            3: 0.92,
            4: 0.94,
            5: 0.95,
            6: 0.96,
            7: 0.965,
            8: 0.97, }

estimator_dilatation_table = defaultdict(lambda: 1)
estimator_dilatation_table.update(pre_dict)


def load_rna_counts_table(rna_source, experiments_to_load):
    """
    Imports the counts tables from a csv file. The csv file has to contain gene identifiers,
    uxon lengths as first two columns and, for each experiment, counts. Lines are genes,
    experiments are columns. In addition to importing, performs conversion to the DB inner IDs
    from heterogeneous identifiers

    :param rna_source: path towards the file where the raw RNA_seq counts are stored
    :param experiments_to_load: number of experiments we want to account retrieve
    """
    gene_names = []
    uxon_lengths = []
    table = []

    print 'loading RNA counts table from %s' % rna_source
    with open(rna_source, 'rb') as source_file:
        rdr = reader(source_file, 'excel-tab')
        rdr.next()  # skipping the headers
        for row in rdr:
            gene_names.append(np.array([row[0], '']))
            counts_for_gene = [float(f) for f in row[2: experiments_to_load + 2]]
            table.append(np.array(counts_for_gene, dtype=np.float64))
            uxon_lengths.append(np.array([float(row[1])]))
    gene_names = np.array(gene_names)
    table = np.array(table)
    uxon_lengths = np.array(uxon_lengths)
    return gene_names, uxon_lengths, table


def counts_filter(table, experiment_groups, filter_level):
    """
    Generates a boolean array that filters out all the groups where counts are below defined level

    :param table: table of counts; line is gene, row is experiment
    :param experiment_groups: groups experiments into repeats
    :param filter_level: minimum amount of counts in any experience for which we are going to
    accept a gene
    :return: a boolean array which is True if Gene passed that test
    """
    minima = np.zeros((table.shape[0], len(experiment_groups)))

    for i, group in enumerate(experiment_groups):
        minima[:, i] = np.min(table[:, group], axis=1)
    filter_mask = np.max(minima, axis=1) > filter_level - 1

    return filter_mask


def convert_to_rpkm(uxon_length, table):
    """
    Translates counts to the RPKMs

    :param uxon_length: vector of gene lengths
    :param table: table of counts
    :return: table of RPKMs
    """
    total_series_reads = np.sum(table, axis=0)
    re_table = table / total_series_reads[np.newaxis, :] / uxon_length[:, 0][:, np.newaxis] * 10e12
    return re_table


def significantly_different_genes(rpkm_table, experiment_groups, intergroups, target_p_value=0.05):
    """
    Performs a test that uses the error function to determine if we can reject the hypothesis that
    all the genes are sampled from the same distribution.

    :param rpkm_table: table of the rpkm values
    :param experiment_groups: groups on indexes
    :param intergroups: the groups between which we want to do the comparisons
    :param target_p_value: p_value with which we want to be able to reject the null hypothesis
    """
    groups_means = np.zeros((rpkm_table.shape[0], len(experiment_groups)))
    groups_var = np.zeros((rpkm_table.shape[0], len(experiment_groups)))

    for i, group in enumerate(experiment_groups):
        groups_means[:, i] = np.mean(rpkm_table[:, group], axis=1)
        groups_var[:, i] = np.var(rpkm_table[:, group], axis=1) / \
            estimator_dilatation_table[len(group)]**2

    group_comparison = []
    for bi_group in intergroups:
        groups_mean_difference = np.fabs(groups_means[:,
                                         bi_group[0]] - groups_means[:, bi_group[1]])
        groups_combined_std = np.sqrt(groups_var[:, bi_group[0]] + groups_var[:, bi_group[1]])
        p_val = t.sf(groups_mean_difference / groups_combined_std, (len(experiment_groups[
            bi_group[0]]) + len(experiment_groups[bi_group[1]]))/2)
        sorted_p_vals = np.sort(p_val, axis=0)
        lower_index = np.array(range(0, sorted_p_vals.shape[0])) *\
            target_p_value / sorted_p_vals.shape[0]
        pre_filter_mask = sorted_p_vals <= lower_index
        filter_mask = pre_filter_mask
        if np.any(pre_filter_mask):
            refined_threshold = np.max(sorted_p_vals[pre_filter_mask])
            filter_mask = p_val < refined_threshold
        group_comparison.append((p_val, filter_mask))

    return group_comparison


def run_analysis_suite(rna_source, no_of_experiments, experimental_groups, groups_to_compare,
                       count_filter_level=5, false_discovery_rate=0.05):
    """
    Imports counts table, runs test suite and stores the result of statistical analysis for further
    computation. returns stored values to the standard output.

    :param rna_source: the file from which the raw counts are to be read
    :param no_of_experiments: number of experiments
    :param experimental_groups: experiment groupings
    :param groups_to_compare: groups to be compared
    :param count_filter_level: minimum counts to run statistics
    :param false_discovery_rate: desired false discovery rate
    """

    names, lengths, counts = load_rna_counts_table(rna_source, no_of_experiments)
    _, _, names[:, 1] = DB_IO_Routines.look_up_annotation_set(names[:, 0].tolist())

    filter_mask = counts_filter(counts, experimental_groups, filter_level=count_filter_level)

    names = names[filter_mask, :]
    lengths = lengths[filter_mask, :]
    counts = counts[filter_mask, :]

    rpkms = convert_to_rpkm(lengths, counts)
    testres = significantly_different_genes(rpkms, experimental_groups,
                                            groups_to_compare, false_discovery_rate)
    filter_masks = [test_[1] for test_ in testres]
    dump_object(Dumps.RNA_seq_counts_compare, filter_masks)

    return filter_masks


if __name__ == "__main__":
    pp = PrettyPrinter(indent=4)

    exp_groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    test_groups_to_compare = [[0, 1], [0, 2]]

    run_analysis_suite(RNA_source, 9, exp_groups, test_groups_to_compare, 10, 0.05)
