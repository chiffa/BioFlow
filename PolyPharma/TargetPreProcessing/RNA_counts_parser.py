__author__ = 'ank'

from PolyPharma.configs import RNA_source, Outputs
import numpy as np
from itertools import combinations as comb
from math import factorial
from csv import reader, writer

# TODO:
    # add filters on the counts (max of mins in classes)
    # add computation of RPKM
    # add computation of average class RPKM
    # add compuatation of RPKM std deviation ESTIMATOR!
    # add computation of intercoverage rate of Gaussians (use 'erf') (i.e. p_value of coming from the same gaussian)


def reshape(data, ms_array):
    """ """
    return np.reshape(ms_array, (data.shape[0],1))


def measure_pass_rate(test_name, result):
    """ """
    print test_name, " passed: ", result.nonzero()[0].__len__()
    return  result.nonzero()[0].__len__()


def list_diff(list_a, list_b):
    """ """
    return set(list_a)-set(list_b)


def import_counts_table(counts_size):
    """
    """
    names = np.zeros((1,1))
    uxon_lenghts = np.zeros((1,1))
    table = np.zeros((counts_size, 1)).T

    with open(RNA_source, 'rb') as source_file:
        rdr = reader(source_file, 'excel-tab')
        rdr.next()
        for row in rdr:
            rrw = [float(f) for f in row[2 : counts_size + 2] ]
            table = np.concatenate( (table, np.array([rrw], dtype = np.float64)))
            uxon_lenghts = np.concatenate( (uxon_lenghts, np.array([row[1: 2]])))
            names = np.concatenate( (names, np.array([row[: 1]])))
    print table.shape
    return table[1: , :], names[1: , :]



def import_data(anot_size, data_size):
    """
    """
    names = np.zeros((anot_size,1), dtype = np.float64).T
    table = np.zeros((data_size, 1)).T
    with open(RNA_source, 'rb') as source_file:
        rdr = reader(source_file, 'excel-tab')
        rdr.next()
        for row in rdr:
            rrw = [float(f) for f in row[anot_size : data_size + anot_size] ]
            table = np.concatenate( (table, np.array([rrw], dtype = np.float64)))
            names = np.concatenate( (names, np.array([row[0 : anot_size]])))
    print table.shape
    return table[1: , :], names[1: , :]


def export_test_results(names, results):
    final_res = names.__copy__()
    for result in results:
        final_res = np.concatenate((final_res, result), axis = 1)
    with open(Outputs.RNA_pre_filter_output, 'wb') as out_file:
        out_writer = writer(out_file, 'excel')
        for i in range(0, final_res.shape[0]):
            out_writer.writerow(final_res[i,:].flatten().tolist())


def prima_filter(data, groups):
    """
    filters out samples where there are too much intra_sample noise
    """
    noise_list = []
    for group in groups:
        std = np.std(data[:, group], axis = 1)
        mean = np.mean(data[:, group], axis = 1)
        noise_list.append(np.reshape(std/mean,(data.shape[0],1)))

    noise = np.concatenate(tuple(noise_list), axis = 1)
    return np.max(noise, axis = 1)


def noise_removed_difference(data, group_1, group_2):
    """ """
    var_1 = np.var(data[:, group_1], axis = 1)
    mean_1 = np.mean(data[:, group_1], axis = 1)
    var_1, mean_1 = (reshape(data, var_1), reshape(data, mean_1))

    var_2 = np.var(data[:, group_2], axis = 1)
    mean_2 = np.mean(data[:, group_2], axis = 1)
    var_2, mean_2 = (reshape(data, var_2), reshape(data, mean_2))

    intra_v = np.mean(np.concatenate((var_1,var_2), axis=1), axis=1)
    inter_v = np.var(np.concatenate((mean_1,mean_2), axis=1), axis=1)

    return reshape(data, inter_v / intra_v), reshape(data, np.log2(mean_2 / mean_1))


def permutation_test(test, data, group_1, group_2):
    """

    """
    # TODO: render it probabilist over some threshold of samples (total indexes > 12 for 1k, > 15 for 10k )
    sh = (data.shape[0], 1)
    total_indexes = group_1 + group_2
    total = factorial(len(total_indexes)) / factorial(len(group_1)) / factorial(len(group_2))
    th_lim = 1/float(total)
    print "theoretical p-val limit for ", str((group_1, group_2)), ": ", th_lim
    counter = np.zeros(sh, dtype = np.float64)
    default_test = test(data, group_1, group_2)[0]
    for g_1 in comb(total_indexes, group_1.__len__()):
        g_2 = tuple(list_diff(total_indexes,g_1))
        counter[test(data, g_1, g_2)[0] <= default_test] += 1.0

    p_val = np.ones(sh, dtype = np.float64) - counter/float(total)
    p_val[p_val < th_lim] = th_lim
    return p_val

def test_suite_1(data, prima_test, group_1, group_2, noise_tolerance, signal_strength, p_val, min_gain):
    """

    """
    p_test = reshape(data, prima_test < noise_tolerance)
    nrd, absg = noise_removed_difference(data, group_1, group_2)
    sig_str = reshape(data, nrd > signal_strength)
    absgain = reshape(data, np.absolute(absg) > min_gain)
    ptrv = permutation_test(noise_removed_difference, data, group_1, group_2)
    per_t_p_val = reshape(data,  ptrv <= p_val)

    final = np.logical_and(absgain,
                np.logical_and(p_test,
                    np.logical_and(sig_str,
                                   per_t_p_val)))

    prime_noise, sig_str, per_t_p_val = ( reshape(data, prima_test),
                                          reshape(data, nrd),
                                          reshape(data, ptrv))

    total = measure_pass_rate("full test_suite", final)

    f_names, f_results = (np.reshape(names[final], (total, 1)),
                          (np.reshape(data[final[:, 0], :], (total, data.shape[1])),
                            np.reshape(prime_noise[final], (total, 1)),
                            np.reshape(sig_str[final], (total, 1)),
                            np.reshape(per_t_p_val[final], (total, 1)),
                            np.reshape(absg[final[:,0]], (total, 1)))
                            )

    return f_names, f_results

if __name__ == "__main__":
    data, names = import_data(1, 9)

    groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    prima_test = prima_filter(data, groups)

    f_names, f_results = test_suite_1(data,prima_test, groups[0], groups[1]+groups[2], 0.1, 2, 0.05, 1.0)
    # 2 in signal_strength corresponds to a probability of 4.2% of mistaking noise for information in our case.
    # However, if the distribution isn't gaussian, it doesn't hold well.

    export_test_results(f_names, f_results)
















# def read_RNA_seq_counts():
#     retlist = []
#     print 'opening: ', RNA_source
#     with open(RNA_source,'rb') as source_file:
#         source_reader = reader(source_file, dialect='excel-tab')
#         for row in source_reader:
#             retlist.append(row)
#     return retlist
#
# def filter_id_only(parse_table):
#     return [sublist[0] for sublist in parse_table if 'ENS' in sublist[0]]
#
#
#
#

#
#
#
#
# if __name__ == "__main__":
#     print filter_id_only(read_RNA_seq_counts())