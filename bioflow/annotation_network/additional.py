from bioflow.annotation_network.knowledge_access_analysis import get_go_interface_instance, \
    ref_param_set, compare_to_blank, log
from bioflow.utils.linalg_routines import analyze_eigenvects
from copy import copy
from random import shuffle


def linearly_independent_go_groups(size):
    """
    Performs the analysis of linearly independent GO groups

    :param size: number of most prominent eigenvectors contributing to an eigenvalue
    we would want to see
    """
    go_interface_instace = get_go_interface_instance(ref_param_set)
    go_interface_instace.undump_independent_linear_sets()
    char_indexes = dict(
        (key,
         (len(
             go_interface_instace.GO2UP_Reachable_nodes[value]),
             go_interface_instace.GO_Legacy_IDs[value],
             go_interface_instace.GO_Names[value])) for key,
        value in go_interface_instace.Num2GO.items())
    print(go_interface_instace.pretty_time())
    analyze_eigenvects(go_interface_instace.Indep_Lapl, size, char_indexes)


def run_analysis(group):
    """
    Performs the whole analysis round retrieving

    :param group:
    :return:
    """
    go_interface_instance = get_go_interface_instance(ref_param_set)
    go_interface_instance.set_uniprot_source(group)
    go_interface_instance.build_extended_conduction_system()
    go_interface_instance.export_conduction_system()
    nr_nodes, nr_groups = compare_to_blank(
        len(group), [1000, 1200], go_interface_instance)
    for group in nr_groups:
        print('run analysis 1: %s' % group)
    for node in nr_nodes:
        print('run analysis 2: %s' % node)


def decide_regeneration():
    """
    A script to decide at what point it is better to recompute a new a network rather
    then go through the time it requires to be upickled.
    The current decision is that for the samples of the size of ~ 100
    reached_uniprots_neo4j_id_list, we are better off unpickling from 4
    and more by factor 2 and by factor 10 from 9
    Previous experiments have shown that memoization with pickling incurred no noticeable
    delay on samples of up to
    50 UPs, but that the storage limit on mongo DB was rapidly exceeded, leading us to
    create an allocated dump file.
    """
    sample_root = []
    rooot_copy = copy(sample_root)
    go_interface_instance = get_go_interface_instance(ref_param_set)
    go_interface_instance.set_uniprot_source(sample_root)
    go_interface_instance.build_extended_conduction_system()
    go_interface_instance.export_conduction_system()
    log.info('decide_regeneration 1: %s', go_interface_instance.pretty_time())
    for i in range(2, 9):
        shuffle(rooot_copy)
        go_interface_instance.export_subsystem(sample_root, rooot_copy[:i ** 2])
        log.info('decide_regeneration 2: %s, retrieve \t %s',
                 i ** 2, go_interface_instance.pretty_time())
        go_interface_instance.set_uniprot_source(rooot_copy[:i ** 2])
        go_interface_instance.build_extended_conduction_system(memoized=False)
        go_interface_instance.export_conduction_system()
        log.info('decide_regeneration 3: %s, redo: \t %s',
                 i ** 2, go_interface_instance.pretty_time())
