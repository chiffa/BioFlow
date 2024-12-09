"""
Stability degradation

"""
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from csv import reader as csv_reader

from collections import defaultdict


def compare_calls(path_1_to_compare, path_2_to_compare,
                  annotome=False, p_val_cutoff=0.0005,
                  experiment_name='', comparison_pair=('', ''),
                  illustrate=True):

    # print(path_1_to_compare, '\n', path_2_to_compare)

    # define the statics

    column_selector = 2  # 2 for both the interactome and annotome for info flow,

    # p_val_cutoff = 0.005
    p_val_column_selector = 4  # 4 for the interactome; 5 for annotome

    if annotome:
        p_val_column_selector = 5

    # group up the flow values in the correct order
    table_1 = np.genfromtxt(path_1_to_compare, delimiter='\t')
    table_2 = np.genfromtxt(path_2_to_compare, delimiter='\t')

    ids_1 = table_1[1:, 0].astype(int)
    ids_2 = table_2[1:, 0].astype(int)

    selection_1 = table_1[1:, column_selector].astype(float)
    selection_2 = table_2[1:, column_selector].astype(float)

    p_val_list_1 = table_1[1:, p_val_column_selector].astype(float)
    p_val_list_1_filter = p_val_list_1 < p_val_cutoff

    p_val_list_2 = table_2[1:, p_val_column_selector].astype(float)
    p_val_list_2_filter = p_val_list_2 < p_val_cutoff

    common_ids = set(ids_1[p_val_list_1_filter].tolist()).intersection(
                 set(ids_2[p_val_list_2_filter].tolist()))

    # finally get the values in the same order from both sets of calls
    accumulate_set_1 = []
    accumulate_set_2 = []

    for id in common_ids:
        accumulate_set_1.append(selection_1[ids_1 == id])
        accumulate_set_2.append(selection_2[ids_2 == id])

    # get spearman correlation
    s_value, s_p_val = spearmanr(np.array(accumulate_set_1),
                                 np.array(accumulate_set_2))

    # print('s/p values: %.2f; %.2f %%' % (s_value, s_p_val*100))

    s_1_excess = len(selection_1[p_val_list_1_filter]) - len(common_ids)
    s_2_excess = len(selection_2[p_val_list_2_filter]) - len(common_ids)

    # print('overlapping, non-overlapping ids in sets 1 & 2: %d, %d & %d' %
          # (len(common_ids), s_1_excess, s_2_excess))

    # now plot the correlation itself.
    plt.title(experiment_name)
    plt.plot(accumulate_set_1, accumulate_set_2, 'ko')
    plt.xlabel(comparison_pair[0])
    plt.ylabel(comparison_pair[1])

    if illustrate:
        plt.show()
    else:
        plt.clf()

    return s_value, s_p_val, len(common_ids), s_1_excess, s_2_excess


def compare_groups(path_1_to_compare, path_2_to_compare,
                   p_val_cutoff=0.0005,
                   experiment_name='', comparison_pair=('', ''),
                   illustrate=True):

    def align_groups(clusters_set):

        # list of the contents of the cluster
        # calculate the angular alignment
        # draw the alignment matrix

        pass

    def parse_groups(path_to_group_files):

        in_group_parse = True
        active_group = -1

        group_contents = defaultdict(list)
        group_chars = {}

        with open(path_to_group_files, 'rt') as source_group:
            reader = csv_reader(source_group, delimiter='\t')
            for line in reader:
                if len(line) < 2:
                    continue
                if line[0] == 'cluster_no':
                    in_group_parse = False
                    continue
                if line[1] == 'id':
                    in_group_parse = True
                    continue
                if not in_group_parse:
                    group_chars[line[0]] = (float(line[1]), int(float(line[2])), float(line[3]))
                    active_group = line[0]
                    continue
                if in_group_parse:
                    group_contents[active_group].append(int(line[1]))

        return group_contents, group_chars

    # print(path_1_to_compare, '\n', path_2_to_compare)

    group_1_contents, group_1_chars = parse_groups(path_1_to_compare)
    group_2_contents, group_2_chars = parse_groups(path_2_to_compare)

    group_1_chars = {key: vals for key, vals in group_1_chars.items() if vals[0] < p_val_cutoff}
    group_2_chars = {key: vals for key, vals in group_2_chars.items() if vals[0] < p_val_cutoff}

    group_1_contents = [(key, vals)
                        for key, vals
                        in group_1_contents .items()
                        if key in group_1_chars.keys()]

    group_2_contents = [(key, vals)
                        for key, vals
                        in group_2_contents .items()
                        if key in group_2_chars.keys()]

    sum_matrix = np.full([len(group_1_contents), len(group_2_contents)], 0.0)

    for i, (key_1, vales_1) in enumerate(group_1_contents):
        for j, (key_2, vales_2) in enumerate(group_2_contents):

            intersect_1 = len(set(vales_1).intersection(set(vales_2))) / len(vales_1)
            intersect_2 = len(set(vales_1).intersection(set(vales_2))) / len(vales_2)

            sum_matrix[i, j] = intersect_1 + intersect_2

    eq_1 = []
    eq_1_p = []
    eq_1_p_vals = []
    eq_1_flow = []
    overflow_1 = []
    overflow_1_flow = []

    eq_2 = []
    eq_2_p = []
    eq_2_p_vals = []
    eq_2_flow =[]
    overflow_2 = []
    overflow_2_flow = []

    if len(group_1_contents) >= len(group_2_contents):

        for i, (key_1, vales_1) in enumerate(group_2_contents):
            eq_2.append(key_1)
            j = np.argmax(sum_matrix[:, i])
            eq_1.append(group_1_contents[int(j)][0])
            eq_2_p.append(i)
            eq_1_p.append(j)

        overflow_1 = list(set(dict(group_1_contents).keys()) - set(eq_1))
        overflow_1_flow = [group_1_chars[key] for key in overflow_1]

    else:

        for i, (key_1, vales_1) in enumerate(group_1_contents):
            eq_1.append(key_1)
            j = np.argmax(sum_matrix[i, :])
            eq_2.append(group_2_contents[int(j)][0])
            eq_1_p.append(i)
            eq_2_p.append(j)

        overflow_2 = list(set(dict(group_2_contents).keys()) - set(eq_2))
        overflow_2_flow = [group_2_chars[key] for key in overflow_2]

    for k_1, k_2 in zip(eq_1, eq_2):

        eq_1_p_vals.append(group_1_chars[k_1][0])
        eq_1_flow.append(group_1_chars[k_1][2])

        eq_2_p_vals.append(group_2_chars[k_2][0])
        eq_2_flow.append(group_2_chars[k_2][2])

        pass

    # print('groups overlap/overflow: %d/%d/%d' % (len(eq_1), len(overflow_1), len(overflow_2)))
    # print('combined flow overlap/overflow: %.2f/%.2f %.2f/%.2f' % (np.sum(np.array(eq_1_flow)),
    #                                                                np.sum(np.array(eq_2_flow)),
    #                                                                np.sum(np.array(overflow_1_flow)),
    #                                                                np.sum(np.array(overflow_2_flow))))

    s_value_pval, s_p_val_pval = spearmanr(np.array(eq_1_p_vals),
                                           np.array(eq_2_p_vals))

    # print('s/p values for p_value: %.2f; %.2f %%' % (s_value_pval, s_p_val_pval*100))

    s_value_flow, s_p_val_flow = spearmanr(np.array(eq_1_flow),
                                           np.array(eq_2_flow))

    # print('s/p values for flow: %.2f; %.2f %%' % (s_value_flow, s_p_val_flow*100))

    plt.title('Cluster alignment matrix for %s' % experiment_name)

    plt.imshow(sum_matrix, interpolation='None', cmap='Reds')
    plt.plot(eq_2_p, eq_1_p, 'kx')

    plt.xlabel(comparison_pair[1])
    plt.ylabel(comparison_pair[0])

    plt.colorbar()

    if illustrate:
        plt.show()
    else:
        plt.clf()

    plt.title('p_val and info flow correlation for %s' % experiment_name)

    ax1 = plt.subplot(1, 2, 2)
    ax1.title.set_text('P value Correlation')
    plt.loglog(eq_2_p_vals, eq_1_p_vals, 'ko')
    plt.xlabel(comparison_pair[0])
    plt.ylabel(comparison_pair[1])

    ax2 = plt.subplot(1, 2, 1)
    ax2.title.set_text('Flow Correlation')
    plt.plot(eq_2_flow, eq_1_flow, 'ko')
    plt.xlabel(comparison_pair[0])
    plt.ylabel(comparison_pair[1])

    if illustrate:
        plt.show()
    else:
        plt.clf()

    return (len(eq_1), len(overflow_1), len(overflow_2)),\
           (s_value_pval, s_p_val_pval),\
           (s_value_flow, s_p_val_flow)



def traversal_comparator(walk_root, anchors_list,
                         reference, illustration=True):

    interactome_call_files = {}
    interactome_group_files = {}

    annotome_call_files = {}
    annotome_group_files = {}

    for dirpath, dirnames, filenames in os.walk(walk_root):
        leaf_folder = os.path.basename(os.path.normpath(dirpath))
        if leaf_folder in anchors_list:
            for file in filenames:
                full_path = os.path.join(dirpath, file)
                if file[-4:] == '.tsv':
                    if 'interactome_analysis' in file:
                        interactome_call_files[leaf_folder] = full_path
                    if 'interactome_clusters' in file:
                        interactome_group_files[leaf_folder] = full_path
                    if 'knowledge_analysis' in file:
                        annotome_call_files[leaf_folder] = full_path
                    if 'knowledge_clusters' in file:
                        annotome_group_files[leaf_folder] = full_path

    for leaf_folder in anchors_list:
        if leaf_folder != reference:
            s_value, s_p_val, len_common_ids, s_1_excess, s_2_excess = \
            compare_calls(interactome_call_files[reference],
                          interactome_call_files[leaf_folder],
                          annotome=False,
                          experiment_name='interactome calls comparison',
                          comparison_pair=(reference, leaf_folder),
                          illustrate=illustration)

            (len_eq_1, len_overflow_1, len_overflow_2),\
            (s_value_pval, s_p_val_pval),\
            (s_value_flow, s_p_val_flow) = \
            compare_groups(interactome_group_files[reference],
                           interactome_group_files[leaf_folder],
                          experiment_name='interactome groups comparison',
                          comparison_pair=(reference, leaf_folder),
                          illustrate=illustration)


            print('Interactome stats for:\t %s\n'
                  '\tCalls:\n'
                  '\t\tCommon/ref overhang/exp overhang:\t %d\t%d\t%d\n'
                  '\t\tCommon: flow spearman_corr/corr p_val: \t %.2f\t%.2f%%\n'
                  '\tGroups:\n'
                  '\t\tCommon/ref overhang/exp overhang:\t %d\t%d\t%d\n'
                  '\t\tCommon: flow spearman_corr/corr p_val: \t %.2f\t%.2f%%\n'
                  '\t\tCommon: p_val spearman_corr/corr p_val: \t %.2f\t%.2f%%\n' %
                  (leaf_folder,
                   len_common_ids, s_1_excess, s_2_excess,
                   s_value, s_p_val*100,
                   len_eq_1, len_overflow_1, len_overflow_2,
                   s_value_flow, s_p_val_flow*100,
                   s_value_pval, s_p_val_pval*100))

            s_value, s_p_val, len_common_ids, s_1_excess, s_2_excess = \
            compare_calls(annotome_call_files[reference],
                          annotome_call_files[leaf_folder],
                          annotome=True,
                          experiment_name='annotome calls comparison',
                          comparison_pair=(reference, leaf_folder),
                          illustrate=illustration)

            (len_eq_1, len_overflow_1, len_overflow_2),\
            (s_value_pval, s_p_val_pval),\
            (s_value_flow, s_p_val_flow) = \
            compare_groups(annotome_group_files[reference],
                           annotome_group_files[leaf_folder],
                          experiment_name='annotome groups comparison',
                          comparison_pair=(reference, leaf_folder),
                          illustrate=illustration)

            print('Annotome stats for:\t %s\n'
                  '\tCalls:\n'
                  '\t\tCommon/ref overhang/exp overhang:\t %d\t%d\t%d\n'
                  '\t\tCommon: flow spearman_corr/corr p_val: \t %.2f\t%.2f%%\n'
                  '\tGroups:\n'
                  '\t\tCommon/ref overhang/exp overhang:\t %d\t%d\t%d\n'
                  '\t\tCommon: flow spearman_corr/corr p_val: \t %.2f\t%.2f%%\n'
                  '\t\tCommon: p_val spearman_corr/corr p_val: \t %.2f\t%.2f%%\n' %
                  (leaf_folder,
                   len_common_ids, s_1_excess, s_2_excess,
                   s_value, s_p_val*100,
                   len_eq_1, len_overflow_1, len_overflow_2,
                   s_value_flow, s_p_val_flow*100,
                   s_value_pval, s_p_val_pval*100))


if __name__ == "__main__":
    path_1_to_compare = 'C:\\Users\\Andrei\\PycharmProjects\\deployments\\BioFlow runs\\' \
                        'run started on 2021-08-10 18.25.32.510454\\' \
                        'kp_km_human\\interactome_analysis_stats.tsv'
    path_2_to_compare = 'C:\\Users\\Andrei\\PycharmProjects\\deployments\\BioFlow runs\\' \
                        'run started on 2021-08-11 23.18.18.729739\\' \
                        'kp_km_human_20p_drop\\interactome_analysis_stats.tsv'

    path_3_to_compare = "C:\\Users\\Andrei\\PycharmProjects\\deployments\\BioFlow runs" \
                        "\\run started on 2021-08-10 18.25.32.510454\\kp_km_human" \
                        "\\interactome_clusters_stats.tsv"
    path_4_to_compare = "C:\\Users\\Andrei\\PycharmProjects\\deployments\\BioFlow runs" \
                        "\\run started on 2021-08-11 23.18.18.729739\\kp_km_human_20p_drop" \
                        "\\interactome_clusters_stats.tsv"

    # compare_calls(path_1_to_compare, path_2_to_compare,
    #               experiment_name="Comparison between humanized kp/km with 20% hits dropped",
    #               comparison_pair=('calls for all hits', 'calls for 20% hits dropped'))

    # compare_groups(path_3_to_compare, path_4_to_compare)

    traversal_comparator('C:\\Users\\Andrei\\PycharmProjects\\deployments\\BioFlow runs\\',
                         ['kp_km_human', 'kp_km_human_20p_drop'],
                         'kp_km_human',
                         illustration=True)