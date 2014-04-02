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
    pass


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
        print sample['sys_hash'], sample['sys_hash'] == MD5_hash
        UP_set = pickle.loads(sample['UPs'])
        Dic_system = KG.compute_and_export_conduction_system(UP_set)
        curr_inf_conf = []
        for value in Dic_system.itervalues():
            if value[1] == 'GO':
                curr_inf_conf.append([float(value[0]), float(value[4]), int(value[5])])
        propmat = np.array(curr_inf_conf)
        print propmat, propmat.shape
        print np.corrcoef(propmat)
        curr_inf_conf_general += curr_inf_conf


if __name__ == "__main__":
    # spawn_sampler(([10, 15], [1, 2]))
    # spawn_sampler_pool(2, [5, 10], [1, 1])
    stats_on_existing_circsys(10)
