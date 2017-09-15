import pickle

import numpy as np
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt

from bioflow.main_configs import interactome_rand_samp
from bioflow.utils.log_behavior import get_logger
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.algorithms_bank.conduction_routines import perform_clustering

log = get_logger(__name__)


interactome_interface_instance = InteractomeInterface(True, False)
interactome_interface_instance.fast_load()

md5_hash = interactome_interface_instance.md5_hash()

print "samples found to test against:\t %s" % interactome_rand_samp.find({'size': 2,
                                                                          'sys_hash': md5_hash,
                                                                          'sparse_rounds': False}).count()

current_accumulator = []
steps_length_accumulator = []

for i, sample in enumerate(interactome_rand_samp.find({'size': 2, 'sys_hash': md5_hash,
                                                                  'sparse_rounds': False})):
    _, nodes_current_dict = pickle.loads(sample['currents'])
    tensions = pickle.loads(sample['voltages'])

    nodes_current = np.sort(np.array(nodes_current_dict.values()).astype(np.float))[-100:]
    node_ids, tension = (tensions.keys()[0], tensions.values()[0])

    # print nodes_current_dict[node_ids[0]], nodes_current_dict[node_ids[1]]

    print '.',
    if nodes_current_dict[node_ids[0]] > 1e-8 and nodes_current_dict[node_ids[1]] > 1e-8:

        print ''
        print 'processing sample #', i
        print node_ids
        print nodes_current_dict[node_ids[0]]
        print nodes_current_dict[node_ids[1]]
        print tension
        print tension/np.mean([nodes_current_dict[node_ids[0]], nodes_current_dict[node_ids[1]]])

        path_length = 1. / (nodes_current / tension)

        nodes_current = nodes_current[path_length > 1]
        path_length = path_length[path_length > 1]

        current_accumulator += nodes_current[path_length > 1].tolist()
        steps_length_accumulator.append(1./tension)



data = np.array(current_accumulator)
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.semilogy(xs, density(xs), 'k')
plt.show()

data = np.array(steps_length_accumulator)
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.plot(xs, density(xs), 'k')
plt.show()


pickle.dump(steps_length_accumulator, open('step_length.dmp', 'w'))