import pickle

import numpy as np
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
from scipy import histogram2d
from csv import reader as csv_reader
import os

from bioflow.main_configs import interactome_rand_samp_db
from bioflow.utils.log_behavior import get_logger
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.algorithms_bank.conduction_routines import perform_clustering
from bioflow.main_configs import Dumps
from bioflow.utils.top_level import map_and_save_gene_ids
# from bioflow.algorithms_bank.conduction_routines import get_current_through_nodes
from matplotlib.cm import get_cmap
from bioflow.main_configs import output_location
from bioflow.user_configs import sources_location


log = get_logger(__name__)


# TODO: for our application a critical validation step could be to cribble the laplacian with zeros
# => That would be an additional validation of our model stability

interactome_interface_instance = InteractomeInterface(True, True)
interactome_interface_instance.fast_load()

md5_hash = interactome_interface_instance.md5_hash()

interactome_interface_instance.randomly_sample(samples_size=[2],
                                               samples_each_size=[1000])

# interactome_interface_instance.export_conduction_system()

# raise Exception('debug')

# we will need to populate the database first
# here as well is where we will be doing filtering by the edge type
print("samples found to test against:\t %s" % interactome_rand_samp_db.count_documents({'size': 2,
                                                                          'sys_hash': md5_hash,
                                                                          'sparse_rounds': False}))

essential_genes_bulbs_ids = []


# underlying file: "C:\\Users\\Andrei\\Dropbox\\workspaces\\EPFL
#       post-doc\\Mehdi_paper_EPFL\\yeast_essential_genes-gene_names.csv"
essential_genes_bulbs_ids, _ = map_and_save_gene_ids(os.path.join(sources_location,
                                                                  'yeast_essential_genes-gene_names.csv'),
        '')

# was needed due to pull from the string formatting
# essential_genes_bulbs_ids = [int(gene) for gene in essential_genes_bulbs_ids]

node_current_values = []
length_width_accumulator = []
essentiality_percentage = []

# we will need to modify the sys_hash to take in account that we are using only one type of
# connections = > Done
for i, sample in enumerate(interactome_rand_samp_db.find({'size': 2,
                                                          'sys_hash': md5_hash,
                                                          'sparse_rounds': False})):

    # if i > 10:
    #     break

    current_acc, nodes_current_dict = pickle.loads(sample['currents'])
    tensions = pickle.loads(sample['voltages'])  # (tension might be broken)
    # print(sorted(nodes_current_dict.items(), key=lambda x: x[1], reverse=True)[:10])

    io_nodes, tension = (list(tensions.items())[0][0], list(tensions.items())[0][1])

    print('debug tension: ', tension)

    # TODO:

    if tension < 0.1:
        continue  # we are hitting near a very tight cluster, so the pathway will be wide and short

    # # Debug section
    # interactome_interface_instance.current_accumulator = current_acc
    # interactome_interface_instance.node_current = nodes_current_dict
    # interactome_interface_instance.UP2UP_voltages = tensions
    # interactome_interface_instance.export_conduction_system()

    # Collects 100 nodes routing most information and cuts the source/sink nodes
    nodes_current = np.sort(
        np.array(list(nodes_current_dict.values())).astype(np.float))[-102:-2] * tension

    # additional filter in case some of the nodes have too small of current
    nodes_current = nodes_current[nodes_current > 0.05]

    # print(nodes_current)

    # not the most efficient implementation, but oh well
    essential_max_current = 0
    for gene in essential_genes_bulbs_ids:
        if nodes_current_dict[gene] / tension > 0.05:
            if nodes_current_dict[gene] > essential_max_current:
                essential_max_current = nodes_current_dict[gene] * tension

    # print(essential_max_current)

    total_resistance = tension  # V=IR and since I=1 for tension calculation, it is resistance
    length_by_width = total_resistance  # R ~ L / W
    print(total_resistance)
    print(total_resistance.dtype)
    print(nodes_current)
    print(nodes_current.dtype)

    shape_characteristic = 1. / nodes_current

    print('\n\n\n>>>>>>>>>>>>>')
    print('sample ', i)
    print('nodes_current', nodes_current)
    print('length/width / tension', length_by_width)
    print('1/mean_width', np.mean(nodes_current))
    # alternative width is max. But in this case we might to remove everything close enough to 0

    # number of branches routing at least 10% of info flow. This is a bit of an oversimplification
    # ideally we need to use an estimator of width - # of major branches routing flow independently
    # Maybe we can use the extreme values distribution based on the bottom X in order to see
    # which one are really the outliers in the network

    if np.any(nodes_current > .1):
        mean_width = 1./np.mean(nodes_current[nodes_current > 0.1])
    else:
        mean_width = 1./np.median(nodes_current)
    length = mean_width * length_by_width

    if length < 1:
        mean_width /= length
        length = 1

    if mean_width < 1:
        length /= mean_width
        mean_width = 1

    if length > 12:
        continue

    print('width', mean_width)
    print('length', length)
    print('essentiality', essential_max_current)
    # print 'io nodes:\t', nodes_current_dict[io_nodes[0]] * tension, nodes_current_dict[io_nodes[1]] * tension

    # print('resistance:\t', total_resistance)
    # print('max current:\t', nodes_current[-1])
    # print('tension:\t', tension)
    # print(nodes_current)

    length_width_accumulator.append((length, mean_width))
    node_current_values += nodes_current.tolist()

    essentiality_percentage.append(min([essential_max_current, 1.]))

    # raise Exception('debug')


    # TODO: debug dump of a couple of pathways into gdf format

    # if any(nodes_current > 1.1):
    #     print nodes_current
    #     raise Exception('debug')

    # print 'tension', tension
    # print 'total_res', total_resistance
    # print 'path_len, nodes_current', np.column_stack((shape_characteristic, nodes_current))

    # if tension > 0.1:   # (if tension is below, we hit a super closely related cluster - 10 direct connections)

        # shape_characteristic = shape_characteristic[shape_characteristic < 20]
        # shape_characteristic = shape_characteristic[shape_characteristic > 0.1]

        # mean_width = 1.
        # mean_length = 1.
        #
        # w_selected = shape_characteristic[shape_characteristic < 1]
        # if np.any(w_selected):
        #     print 'estimated width distribution:\t', np.mean(1. / w_selected)
        #     mean_width = np.mean(1. / w_selected)
        #     width_accumulator += (1. / w_selected).tolist()
        #     if mean_width < 1:
        #         raise Exception('unexpected width mean')
        # else:
        #     print 'estimated width distribution:\t', 1.
        #
        # l_selected = shape_characteristic[shape_characteristic > 1]
        # if np.any(l_selected):
        #     print 'estimated length distribution:\t', np.mean(l_selected)
        #     mean_length = np.mean(l_selected)
        #     length_accumulator += l_selected.tolist()
        # else:
        #     print 'estimated length distribution: \t', 1.
        #
        # # print essential_max_current
        # # print nodes_current[-1]
        # print "essentiality percentage :\t", essential_max_current/nodes_current[-1]*tension
        # essentiality_percentage.append(essential_max_current/nodes_current[-1]*tension)
        # if essential_max_current*tension > nodes_current[-1]:
        #     print essential_max_current
        #     print nodes_current
        #
        # length_width_accumulator.append([mean_length, mean_width])



node_current_values = np.array(node_current_values)
data = node_current_values[node_current_values > 0.0002]
print(data)
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.plot(xs, density(xs), 'k')
plt.xlabel('pathway shape parameter')
plt.ylabel('density of distribution')
# plt.show()
plt.savefig(os.path.join(output_location, 'pathway_shape_parameter.png'))
plt.clf()


_length = np.array(length_width_accumulator)[:, 0]

print('average length', np.mean(_length[_length > 1.99]))
print('std length', np.std(_length[_length > 1.99]))

data = np.array(_length[_length > 1.99])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Length distribution of non-trivial pathways')
plt.plot(xs, density(xs), 'k')
plt.xlabel('length of the pathway')
plt.ylabel('density of distribution')
# plt.show()
plt.savefig(os.path.join(output_location, 'pathway_length_distribution.png'))
plt.clf()


_width = np.array(length_width_accumulator)[:, 1]

print('average width', np.mean(_width[_length > 1.99]))
print('std width', np.std(_width[_length > 1.99]))

data = np.array(_width[_length > 1.99])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Width distribution of non-trivial pathways')
plt.plot(xs, density(xs), 'k')
plt.xlabel('width of the pathway')
plt.ylabel('density of distribution')
plt.savefig(os.path.join(output_location, 'density_estimation_distribution.png'))
plt.clf()
# plt.show()


data = np.array(np.array(essentiality_percentage)[_length > 1.1])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Percentage of pathway throughput lost in case of essential gene deletion')
plt.plot(xs, density(xs), 'k')
plt.axvline(0.7)
plt.xlabel('percentage of current routed through essential genes')
plt.ylabel('density of distribution')
plt.savefig(os.path.join(output_location, 'essential_percentage_distribution.png'))
plt.clf()
# plt.show()


def better2D_desisty_plot(xdat, ydat, thresh=3, bins=(100, 100)):
    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    distortion = (xyrange[1][1] - xyrange[1][0]) / \
        (xyrange[0][1] - xyrange[0][0])
    xdat = xdat * distortion

    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    hh, locx, locy = histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    # values of the histogram where the points are
    hhsub = hh[posx[ind] - 1, posy[ind] - 1]
    xdat1 = xdat[ind][hhsub < thresh]  # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan  # fill the areas with low density by NaNs

    plt.imshow(
        np.flipud(
            hh.T),
        cmap='jet',
        extent=np.array(xyrange).flatten(),
        interpolation='none')
    plt.plot(xdat1, ydat1, '.')


# better2D_desisty_plot(np.array(length_width_accumulator)[:, 0],
#                       np.array(length_width_accumulator)[:, 1])

cm = get_cmap('seismic')

plt.scatter(np.array(length_width_accumulator)[:, 0], np.array(length_width_accumulator)[:, 1],
            c=np.array(essentiality_percentage), cmap=cm, vmin=0., vmax=1.)
plt.colorbar(label='essentiality')
plt.axvspan(0, 2, facecolor='0.5', alpha=0.5)
plt.xlabel('length of the pathway')
plt.ylabel('width of the pathway')
plt.show()
plt.savefig(os.path.join(output_location, 'essentiality_length_vs_width.png'))
plt.clf()

# pickle.dump(length_accumulator, open('step_length.dmp', 'wb'))
# pickle.dump(width_accumulator, open('width_length.dmp', 'wb'))
w_l_accumulator_location = os.path.join(output_location, 'BOWN_NN_w_l_accumulatror.dmp')
pickle.dump(length_width_accumulator, open(w_l_accumulator_location, 'wb'))