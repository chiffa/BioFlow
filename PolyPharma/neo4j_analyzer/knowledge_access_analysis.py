"""
Builds on the knowledge
"""
__author__ = 'ank'

from Matrix_Interactome_DB_interface import  MatrixGetter
from Matrix_Knowledge_DB_Interface import GO_Interface
from multiprocessing import Pool
from PolyPharma.configs import UP_rand_samp, UP_store
from pprint import PrettyPrinter
import pickle
import numpy as np
import matplotlib.pyplot as plt


pprinter = PrettyPrinter(indent=4)
MG = MatrixGetter(True, False)
MG.fast_load()


def spawn_sampler(sample_size_list_plus_iteration_list):
    """
    Spawns a sampler initalized from the default GO_Interface

    :param sample_size_list_plus_iteration_list: combined list of sample swizes and iterations (requried for Pool.map usage)
    """
    print sample_size_list_plus_iteration_list
    filtr = ['biological_process']
    KG = GO_Interface(filtr, MG.Uniprots, 1, True, 3)
    KG.load()
    print KG.pretty_time()
    sample_size_list = sample_size_list_plus_iteration_list[0]
    iteration_list = sample_size_list_plus_iteration_list[1]
    KG.randomly_sample(sample_size_list, iteration_list)


def spawn_sampler_pool(pool_size, sample_size_list, interation_list_per_pool):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interation_list_per_pool: number of iterations performing the pooling of the samples in each list
    """
    p = Pool(pool_size)
    payload = [(sample_size_list, interation_list_per_pool)]
    p.map(spawn_sampler, payload * pool_size)


def show_corrs(tri_corr_array):
    plt.figure()

    plt.subplot(331)
    plt.title('closest side to side')
    plt.hist(tri_corr_array[0, :], bins=100, histtype='step')

    plt.subplot(332)
    plt.title('closest side to side')
    plt.scatter(tri_corr_array[0, :], tri_corr_array[1, :])

    plt.subplot(333)
    plt.title('closest side to side')
    plt.scatter(tri_corr_array[0, :], tri_corr_array[2, :])

    plt.subplot(335)
    plt.title('closest side to side')
    plt.hist(tri_corr_array[1, :], bins=100, histtype='step')

    plt.subplot(336)
    plt.title('closest side to side')
    plt.scatter(tri_corr_array[1, :], tri_corr_array[2, :])

    plt.subplot(339)
    plt.title('closest side to side')
    plt.hist(tri_corr_array[2, :], bins=100, histtype='step')

    plt.show()

    return np.corrcoef(tri_corr_array)


def stats_on_existing_circsys(size):
    """
    Recovers the statistics on the existing circulation systems.

    :return:
    """
    filtr = ['biological_process']
    KG = GO_Interface(filtr, MG.Uniprots, 1, True, 3)
    KG.load()
    MD5_hash = KG._MD5hash()

    curr_inf_conf_general = []
    for sample in UP_rand_samp.find({'size':size, 'sys_hash' : MD5_hash}):
        UP_set = pickle.loads(sample['UPs'])
        curr_mat, node_currs = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        Dic_system = KG.compute_conduction_system(curr_mat,UP_set,node_currs)
        curr_inf_conf = []
        propmat = np.array(list(Dic_system.itervalues())).T
        curr_inf_conf_general += curr_inf_conf
    print show_corrs(curr_inf_conf_general)


if __name__ == "__main__":
    # spawn_sampler(([10, 15], [1, 2]))
    spawn_sampler_pool(6, [5, 10], [10, 10])
    # stats_on_existing_circsys(10)
