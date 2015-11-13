"""
New analytical routines for the interactome
"""
__author__ = 'ank'
import pickle
import numpy as np
from collections import namedtuple
from multiprocessing import Pool
from pprint import PrettyPrinter
from matplotlib import pyplot as plt
from csv import reader
from Matrix_Interactome_DB_interface import  MatrixGetter
from BioFlow.configs2 import Interactome_rand_samp, Dumps, prename2
from BioFlow.Utils.dataviz import kde_compute
from BioFlow.neo4j_analyzer.Conduction_routines import perform_clustering
from BioFlow.neo4j_analyzer.IO_Routines import undump_object


def MG_gen():
    """
    Retrieves an "interactome_interface" object

    :return:
    """
    MG = MatrixGetter(True, False)
    MG.fast_load()
    print MG.pretty_time()
    return MG


def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initalized from the default GO_Interface

    :param sample_size_list_plus_iteration_list: combined list of sample swizes and iterations (requried for Pool.map usage)
    """
    MG_object = sample_size_list_plus_iteration_list_plus_args[4]
    if MG_object is None:
        MG = MG_gen()
    else:
        MG = MG_object

    sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
    iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
    sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
    chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]
    MG.randomly_sample(sample_size_list, iteration_list, sparse_rounds, chromosome_specific)


def spawn_sampler_pool(pool_size, sample_size_list, interation_list_per_pool, sparse_rounds=False, chromosome_specific=False, MG_object=None):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interation_list_per_pool: number of iterations performing the pooling of the samples in each list
    :param sparse_rounds:
    :type sparse_rounds: int
    :param chromosome_specific:
    :type chromosome_specific: int
    """
    p = Pool(pool_size)
    payload = [(sample_size_list, interation_list_per_pool, sparse_rounds, chromosome_specific, MG_object)]
    p.map(spawn_sampler, payload * pool_size)


def select(bi_array, array_column, selection_span):
    """
    Convinient small function to select a from tri_array all the elements where the column number array_column
    is within the selection span

    :param bi_array: the matrix on which we will be performing the selection
    :param array_column: column on which the selection span will be applied
    :param selection_span: span for which we are going to keep the column.
    :return:
    """
    selector = np.logical_and(selection_span[0]< bi_array[array_column, :], bi_array[array_column, :]< selection_span[1])
    if not any(selector):
        return np.array([[0.0,0.0,0.0]])
    decvec = bi_array[:, selector]
    return decvec


def show_corrs(bi_corr_array, meancorrs, eigvals, selector, test_bi_corr_array, test_meancorr, eigval, resamples):

    # TODO: there is a lot of repetition depending on which values are the biggest, test-setted or real setted.
    # In all, we should be able to reduce it to two functions: scatterplot and histogram with two sets that
    # should go into the dataviz module

    """
    A general function that performs demonstration of an example of random samples of the same size as our sample
    and of our sample and conducts the statistical tests on wherther any of nodes or functional groups in our
    sample are non-random

    :param bi_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the random samples
    :param meancorrs: [[cluster size, average internode connection], ...] - characteristics of clustering random samples with the same parameters
    :param eigvals: eigenvalues associated to the interconnection matrix of random samples
    :param selector: range on which we would like to visually zoom and plot a histogram
    :param test_bi_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the true sample. If none, nothing happens
    :param test_meancorr: [[cluster size, average internode connection], ...] - characteristics of clustering the true sample
    :param eigval: eigenvalues associated to the interconnection matrix of the true sample
    :param resamples: how many random samples we analyzed for the default model
    :return:
    """
    plt.figure()

    plt.subplot(331)
    plt.title('current through nodes')
    bins = np.linspace(bi_corr_array[0, :].min(), bi_corr_array[0, :].max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[0, :].min(), test_bi_corr_array[0, :].min()),
                           max(bi_corr_array[0, :].max(), test_bi_corr_array[0, :].max()),
                           100)
    plt.hist(bi_corr_array[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[0, :], bins=bins, histtype='step', log=True, color='r')

    plt.subplot(332)
    plt.title('test current vs degree')
    plt.scatter(bi_corr_array[1, :], bi_corr_array[0, :])
    if test_bi_corr_array is not None:
        plt.scatter(test_bi_corr_array[1, :], test_bi_corr_array[0, :], color='r' , alpha=0.5)
    plt.axvspan( selector[0], selector[1], facecolor='0.5', alpha=0.3)

    plt.subplot(333)
    plt.title('Currently empty')

    plt.subplot(334)
    plt.title('Gaussian KDE current_info')
    estimator_function = kde_compute(bi_corr_array[(1,0), :], 50, resamples)
    current_info_rel = None
    if test_bi_corr_array is not None:
        current_info_rel = estimator_function(test_bi_corr_array[(1,0),:])

    plt.subplot(335)
    plt.title('Node degree distribution')
    bins = np.linspace(bi_corr_array[1, :].min(), bi_corr_array[1, :].max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[1, :].min(), test_bi_corr_array[1, :].min()),
                           max(bi_corr_array[1, :].max(), test_bi_corr_array[1, :].max()),
                           100)
    plt.hist(bi_corr_array[1, :], bins=100, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[1, :], bins=100, histtype='step', log=True, color='r')

    plt.subplot(336)
    plt.title('Density of current in the highlighted area')
    if test_bi_corr_array is not None:
        bins = np.linspace(min(select(bi_corr_array, 1, selector)[0, :].min(),
                               select(test_bi_corr_array, 1, selector)[0, :].min()),
                           max(select(bi_corr_array, 1, selector)[0, :].max(),
                               select(test_bi_corr_array, 1, selector)[0, :].max()),
                           100)
    plt.hist(select(bi_corr_array, 1, selector)[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(select(test_bi_corr_array, 1, selector)[0, :], bins=100, histtype='step', log=True, color='r')

    # this property is better off viewed as a scatterplot of true points and default points
    plt.subplot(337)
    plt.title('Clustering correlation')
    # plt.scatter(meancorrs[0, :], meancorrs[1, :], color = 'b')
    estimator_function = kde_compute(meancorrs[(0,1), :], 50, resamples)
    cluster_props =None
    if test_meancorr is not None:
        plt.scatter(test_meancorr[0,:], test_meancorr[1,:], color = 'k', alpha=0.8)
        cluster_props = estimator_function(test_meancorr[(0,1),:])

    plt.subplot(338)
    plt.title('Eigvals_hist')
    bins = np.linspace(eigvals.min(), eigvals.max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(eigvals.min(), eigval.min()),
                           max(eigvals.max(), eigval.max()),
                           100)
    plt.hist(eigvals, bins=bins, histtype='step', color = 'b')
    if eigval is not None:
        plt.hist(eigval.tolist()*3, bins=bins, histtype='step', color = 'r')

    plt.subplot(339)
    plt.title('Currently empty')


    plt.show()

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


def compare_to_blanc(blanc_model_size, zoom_range_selector, real_interactome_interface = None, p_val=0.05, sparse_rounds=False, clusters=3, MG_Object=None):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blanc_model_size: the number of uniprots in the blanc model
    :param zoom_range_selector: tuple representing the coverage range for which we would want to see the histogram of current distributions
    :param real_interactome_interface: The interactome_Interface that has run the current computation
    :param p_val: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used with the number of rounds equal to the number
    :type sparse_rounds: int
    :param clusters: specifies the number of clusters we want to have
    :return: None if no significant nodes, the node and group characterisitc dictionaries otherwise
    """
    if MG_Object is None:
        MG = MatrixGetter(True, False)
        MG.fast_load()
    else:
        MG = MG_Object

    MD5_hash = MG._MD5hash()

    curr_inf_conf_general = []
    count = 0
    meancorr_acccumulator = []
    eigval_accumulator = []

    print "samples found to test against:\t", Interactome_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds}).count()
    # this part computes the items required for the creation of a blanc model
    for i, sample in enumerate(Interactome_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds})):
        if sparse_rounds:
            raise Exception('blanc done on sparse rounds, clustering likely to be hazardous. Process interrupted')
        _, node_currs = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        _, _, meancorr, eigvals = perform_clustering(tensions, clusters, show=False)
        meancorr_acccumulator.append(np.array(meancorr))
        eigval_accumulator.append(eigvals)
        Dic_system = MG.format_Node_props(node_currs)
        curr_inf_conf = list(Dic_system.itervalues())
        curr_inf_conf_general.append(np.array(curr_inf_conf).T)
        count = i

    # This part declares the pre-operators required for the verification of a real sample

    final = np.concatenate(tuple(curr_inf_conf_general), axis=1)
    final_meancorrs = np.concatenate(tuple(meancorr_acccumulator), axis=0).T
    final_eigvals = np.concatenate(tuple(eigval_accumulator), axis=0).T
    curr_inf_conf = None
    meancorr = None
    eigval = None
    group2avg_offdiag = None
    node_ids = None
    Dic_system = None
    if real_interactome_interface:
        node_currs = real_interactome_interface.node_current
        Dic_system = MG.format_Node_props(node_currs)
        curr_inf_conf_tot = np.array([[int(key)] + list(val) for key, val in Dic_system.iteritems()]).T
        node_ids, curr_inf_conf = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1,2), :])
        group2avg_offdiag, _, meancorr, eigval = perform_clustering(real_interactome_interface.UP2UP_voltages, clusters)

    print "stats on %s samples" % count


    # TODO: We could and should separate the visualisation from the gaussian estimators computation
    r_nodes, r_groups = show_corrs(final, final_meancorrs, final_eigvals, zoom_range_selector, curr_inf_conf, meancorr.T, eigval.T, count)


    Interactome_node_char = namedtuple('Node_Char',['name', 'current', 'degree', 'p_value'])
    Group_char = namedtuple('Group_Char', ['UPs','num_UPs','average_connection','p_value'])

    if r_nodes is not None:
        not_random_nodes = [str(int(node_id)) for node_id in node_ids[r_nodes < p_val].tolist()]
        not_random_groups = np.concatenate((group2avg_offdiag, np.reshape(r_groups, (3,1))), axis=1)[r_groups < p_val].tolist()
        not_random_groups = [ Group_char(*nr_group) for nr_group in not_random_groups]
        # basically the second element below are the nodes that contribute to the information flow through the node that is considered as
        # non-random
        dct = dict((nr_node_id, Interactome_node_char(MG.ID2displayName[nr_node_id],*(Dic_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist())))
                     for nr_node_id in not_random_nodes)

        return  sorted(dct.iteritems(), key=lambda x:x[1][3]), not_random_groups

    return  None, None, None

# TODO: stabilize the background behavior with regard to the buffering
def auto_analyze(sourcelist, depth, processors=4, backgroundlist=None):
    """
    Automatically analyzes the itneractome synergetic action of the RNA_seq results

    """
    # noinspection PyTypeChecker
    for lst in sourcelist:
        print lst, len(lst)
        MG1 = MG_gen()
        MG1.set_Uniprot_source(list(lst))
        MG1.background = backgroundlist
        print len(MG1.analytic_Uniprots)
        if len(MG1.analytic_Uniprots) < 200:
            spawn_sampler_pool(processors, [len(MG1.analytic_Uniprots)], [depth], MG_object=MG1)
            MG1.build_extended_conduction_system()
            nr_nodes, nr_groups = compare_to_blanc(len(MG1.analytic_Uniprots), [0.5, 0.6], MG1, p_val=0.9, MG_Object=MG1)
        else:
            sampling_depth = max(200**2/len(MG1.analytic_Uniprots), 5)
            print 'lenght: %s \t sampling depth: %s \t, estimated_time: %s' % (len(MG1.analytic_Uniprots), sampling_depth, len(MG1.analytic_Uniprots)*sampling_depth/2/depth/60)
            spawn_sampler_pool(processors, [len(MG1.analytic_Uniprots)], [depth], sparse_rounds=sampling_depth, MG_object=MG1)
            MG1.build_extended_conduction_system(sparse_samples=sampling_depth)
            MG1.export_conduction_system()
            nr_nodes, nr_groups = compare_to_blanc(len(MG1.analytic_Uniprots), [0.5, 0.6], MG1, p_val=0.9, sparse_rounds=sampling_depth, MG_Object=MG1)

        MG1.export_conduction_system()
        for group in nr_groups:
            print group
        for node in nr_nodes:
            print node


def get_source():
    retlist=[]
    with open(prename2) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row
    retlist = [ret for ret in retlist]
    return retlist


if __name__ == "__main__":

    pprinter = PrettyPrinter(indent=4)
    MG = MatrixGetter(True, False)
    MG.fast_load()

    # dumplist = undump_object(Dumps.RNA_seq_counts_compare)

    transcription = ['1005777', '1001874', '842142', '836106', '1014143', '1058552', '821021', '826066', '886586',
                     '865200', '835714', '766912', '900540', '811360', '816278', '1079377', '740496', '1002986',
                     '778615', '1000616', '950906', '996907', '1033828', '971819', '809996', '781060', '874882', '758558',
                     '740148', '758263', '926744', '907135', '990196', '743186', '1011403', '979456', '1073345', '913095',
                     '805935', '777825', '1028447', '951846', '1075746', '882269', '888669', '1029555', '941071', '751245',
                     '986669', '927801', '783375', '1066551', '797230', '859473', '914620', '1083041', '837861', '901826',
                     '906807', '913721', '991031', '831309', '810280', '856443', '986631', '800198', '832809', '774804',
                     '1048252', '935385', '920480', '952510', '988829', '744318', '1001042', '1069894', '848924', '949209',
                     '1035740', '770322', '749382', '1003600', '889200', '1071170', '841076', '822788', '810700', '801219',
                     '801191', '964332', '935139', '1071008', '1002847', '1022239', '875186', '848997', '966988', '1008725',
                     '899845', '763208', '823114', '1078022', '969713', '892913', '1104135', '794981', '1005373', '767473',
                     '820454', '1066442', '753200', '1000686', '800665', '991315', '817329', '885986', '1019819', '871354',
                     '1073014', '886044', '1004612', '871049', '1034170', '860442', '848198', '878608', '844231', '794289',
                     '1016545', '1023218', '1103896', '785577', '995803', '912780', '897083', '886840', '750127', '945831',
                     '916664', '864667', '746441', '1039958', '1101959', '956599', '875904', '1073318', '909015', '974199',
                     '1037879', '964653', '784732', '949392', '1091093', '1081714', '837765', '973051', '772194', '765760',
                     '973374', '893485', '838733', '942660', '798235', '996784', '914154', '1055665', '935238', '780699',
                     '1016715', '991271', '955633', '895017', '794788', '1091063', '1046353', '986976', '826048', '859213',
                     '839038', '783622', '760894', '853760', '1106249', '819917', '1062291', '815494', '862693', '814921',
                     '938485', '783138', '885381', '903943', '772264', '980817', '823252', '962633', '758463', '755526',
                     '757978', '765208', '863431', '1042302', '905303', '811490', '939476', '846291', '993567', '819257',
                     '886290', '1013510', '827509', '991517', '1026799', '766799', '1029295', '806975', '911764', '1070622',
                     '891117', '782645', '994623', '1013051', '830556', '844362', '1040267', '933951', '876062', '1060754',
                     '826120', '973468', '949926', '789918', '747228', '764653', '1056009', '893389', '828460', '780732',
                     '1061892', '748295', '890461', '800707', '896612', '1061111', '1038976', '1026861', '892596', '994512',
                     '908624', '819172', '949891', '1015019', '802349', '752772', '1041116', '920885', '851881', '758797',
                     '898171', '844767', '775297', '936117', '1089489', '970498', '854601', '932389', '911677', '773272',
                     '964086', '837996', '968860', '954895', '847543']


    translation = ['860713', '1054645', '1050811', '1010442', '746034', '960074', '841626', '989237', '992333', '885431',
                   '1049088', '1005233', '742678', '1025641', '1081901', '940598', '1004992', '738650', '1018277',
                   '997379', '745984', '1051693', '758243', '1085622', '933193', '1070658', '932056', '774142',
                   '1001023', '1030061', '743020', '905602', '998094', '932332', '996254', '823270', '856693',
                   '1017722', '976875', '960939', '890389', '1099836', '829293', '957345', '1020969', '905383',
                   '1009002', '771146', '1081329', '985624', '904324', '1019095', '937024', '1051946', '854456',
                   '1076560', '741935', '1102657', '1021881', '1041371', '757039', '997830', '942299', '773580',
                   '902901', '1069785', '1023488', '1077664', '864199', '815792', '947643', '983830', '1049137',
                   '821120', '1094403', '958096', '877347', '870065', '779276', '998107', '966557', '876290', '782349',
                   '886637', '973293', '943484', '840896']


    # MG1.randomly_sample([150], [1], chromosome_specific=15, No_add=True)
    # nr_nodes, nr_groups = compare_to_blanc(150, [0.5, 0.6], MG1, p_val=0.9)
    # MG1.export_conduction_system()
    # for group in nr_groups:
    #     print group
    # for node in nr_nodes:
    #     print node

    #TODO: add loading method for the analysis of the interactome

    sourc = get_source()
    # print sourc

    auto_analyze([sourc], 20)
    pass