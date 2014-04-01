"""
Builds on the knowledge
"""
__author__ = 'ank'

from Matrix_Interactome_DB_interface import  MatrixGetter
from Matrix_Knowledge_DB_Interface import GO_Interface
from multiprocessing import Pool
from PolyPharma.configs import UP_rand_samp, UP_store

MG = MatrixGetter(True, False)
MG.fast_load()



def spawn_sampler(sample_size_list_plus_iteration_list):
    """
    Spawns a sampler initalized from the default GO_Interface

    :param sample_size_list_plus_iteration_list: combined list of sample swizes and iterations (requried for Pool.map usage)
    """
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



def stats_on_existing_circsys(size):
    """

    :return:
    """
    filtr = ['biological_process']
    KG = GO_Interface(filtr, MG.Uniprots, 1, True, 3)
    KG.load()
    MD5_hash = KG._MD5hash()

    for sample in UP_rand_samp.find({'size':size, 'hash':MD5_hash}):
        UP_set = sample['UPs']
        KG.compute_and_export_conduction_system(UP_set)
        # TODO: finish the matrix of knowledge and then the method itself


if __name__ == "__main__":
    spawn_sampler_pool(4, [100, 150, 200, 500], [5, 5, 2, 1])

