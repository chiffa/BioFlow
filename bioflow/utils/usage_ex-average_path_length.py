import pickle

import numpy as np
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
from scipy import histogram2d
from csv import reader as csv_reader

from bioflow.main_configs import interactome_rand_samp_db
from bioflow.utils.log_behavior import get_logger
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.algorithms_bank.conduction_routines import perform_clustering
from bioflow.main_configs import Dumps
# from bioflow.algorithms_bank.conduction_routines import get_current_through_nodes
from matplotlib.cm import get_cmap


log = get_logger(__name__)


interactome_interface_instance = InteractomeInterface(True, True)
interactome_interface_instance.fast_load()

md5_hash = interactome_interface_instance.md5_hash()

print "samples found to test against:\t %s" % interactome_rand_samp_db.find({'size': 2,
                                                                          'sys_hash': md5_hash,
                                                                          'sparse_rounds': False}).count()

essential_genes_bulbs_ids = []

with open(Dumps.analysis_set_bulbs_ids, 'r') as source:
    reader = csv_reader(source)
    for line in reader:
        essential_genes_bulbs_ids += line

essential_genes_bulbs_ids = [int(gene) for gene in essential_genes_bulbs_ids]

values = []
length_width_accumulator = []
essentiality_percentage = []

for i, sample in enumerate(interactome_rand_samp_db.find({'size': 2, 'sys_hash': md5_hash,
                                                                  'sparse_rounds': False})):

    # if i > 10:
    #     break

    _, nodes_current_dict = pickle.loads(sample['currents'])
    tensions = pickle.loads(sample['voltages'])

    io_nodes, tension = (tensions.keys()[0], tensions.values()[0])
    # this actually should be a multiplication - we divide to normalize to 1 volt, after counting for 1 amp

    nodes_current = np.sort(np.array(nodes_current_dict.values()).astype(np.float))[-100:] * tension

    # not the most efficient implementation, but oh well
    essential_max_current = 0
    for gene in essential_genes_bulbs_ids:
        if nodes_current_dict[gene]/tension > 0.05:
            if nodes_current_dict[gene] > essential_max_current:
                essential_max_current = nodes_current_dict[gene] * tension

    # delete nodes close to 1 (IO)
    # ivide by 2

    if tension > .2:

        total_resistance = tension  # yeah, tension is just a resistance in this context - my labeling error
        length_by_width = total_resistance
        # nodes_current = nodes_current[np.logical_not(np.isclose(nodes_current, np.ones(nodes_current.shape), rtol=1e-03))]/2
        nodes_current = nodes_current[nodes_current < 0.999]
        shape_characteristic = 1. / nodes_current

        print '\n\n\n>>>>>>>>>>>>>'
        print 'sample ', i
        print 'length/width', length_by_width
        # alternative width is max. But in this case we might to remove everything close enough to 0
        # mean_width = 1./np.mean(nodes_current[nodes_current > 0.1])
        mean_width = 1. / np.mean(nodes_current[nodes_current > 0.2])
        length = mean_width * length_by_width

        if length < 1:
            mean_width /= length
            length = 1

        if mean_width < 1:
            length /= mean_width
            mean_width = 1

        print 'width', mean_width
        print 'length', length
        print 'essentiality', essential_max_current
        # print 'io nodes:\t', nodes_current_dict[io_nodes[0]] * tension, nodes_current_dict[io_nodes[1]] * tension

        # print 'resistance:\t', total_resistance
        # print 'max current:\t', nodes_current[-1]
        # print 'tension:\t', tension
        # print nodes_current

        length_width_accumulator.append((length, mean_width))
        values += nodes_current.tolist()

        essentiality_percentage.append(min([essential_max_current,1.]))


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


values = np.array(values)

data = values[values > 0.2]
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.plot(xs, density(xs), 'k')
plt.xlabel('pathway shape parameter')
plt.ylabel('density of distribution')
plt.show()


_length = np.array(length_width_accumulator)[:, 0]

print 'average length', np.mean(_length[_length > 1.99])
print 'std length', np.std(_length[_length > 1.99])

data = np.array(_length[_length > 1.99])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Length distribution of non-trivial pathways')
plt.plot(xs, density(xs), 'k')
plt.xlabel('length of the pathway')
plt.ylabel('density of distribution')
plt.show()

_width = np.array(length_width_accumulator)[:, 1]

print 'average width', np.mean(_width[_length > 1.99])
print 'std width', np.std(_width[_length > 1.99])

data = np.array(_width[_length > 1.99])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Width distribution of non-trivial pathways')
plt.plot(xs, density(xs), 'k')
plt.xlabel('width of the pathway')
plt.ylabel('density of distribution')
plt.show()

data = np.array(np.array(essentiality_percentage)[_length > 1.1])
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Percentage of pathway throughput lost in case of essential gene deletion')
plt.plot(xs, density(xs), 'k')
plt.axvline(0.7)
plt.xlabel('percentage of current routed through essential genes')
plt.ylabel('density of distribution')
plt.show()

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

# pickle.dump(length_accumulator, open('step_length.dmp', 'w'))
# pickle.dump(width_accumulator, open('width_length.dmp', 'w'))
pickle.dump(length_width_accumulator, open('w_l_accumulator.dmp', 'w'))