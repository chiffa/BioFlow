import numpy as np
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch


def compute_tension_clustering(voltage_pair_dict: dict,
                               random_sample: bool = True) -> (list, np.array, np.array):
    """
    When supplied with the potential differences, clusters hierarchically to get groups of nodes
    with the highest information flow in between them and returns the list of nodes per cluster
    number, lowest internode information flow in the cluster and number of nodes in the cluster.

    :param voltage_pair_dict: {(internal_id, internal_id): potential diff} basically UP2UP_voltages
    :param random_sample: if true, cluster membership will not be computed.
    :return:what each cluster contains ([] if random_sample is True), average information flow
        between nodes and nodes in cluster
    """

    ups = set()

    for up_1, up_2 in voltage_pair_dict.keys():
        ups.add(up_1)
        ups.add(up_2)

    ups = list(ups)
    total_nodes = len(ups)

    distmat = np.zeros((total_nodes, total_nodes))

    id_2_idx = {up: i for i, up in enumerate(ups)}
    idx_2_id = {i: up for i, up in enumerate(ups)}

    for (up_1, up_2), value in voltage_pair_dict.items():
        idx_1, idx_2 = (id_2_idx[up_1], id_2_idx[up_2])
        distmat[idx_1, idx_2] = value
        distmat[idx_2, idx_1] = value

    cond_distmat = ssd.squareform(distmat)
    clust_linkmap = sch.linkage(cond_distmat)

    # clust_linkmap is basically:
    # [[clust_no_1, clust_no_2, clust_dist, nodes_in_clust]]
    # where clust_no < nodes_no are individual nodes
    # and each new merger i generates a cluster with index n+i, that will be used in further merges
    #
    #
    # what we want is to calculate the map of ids and node nums to internode flows
    # and have a list of nodes, so that if its statistically significant, we can just discard it
    # the simplest thing would be to simply strip the last two columns for the test sets,

    if random_sample:
        # [], average information flow between nodes and nodes in cluster
        return [], 1./clust_linkmap[:, 2], clust_linkmap[:, 3]

    else:
        clusters_collection = []

        for line in clust_linkmap:
            id_1, id_2 = (line[0], line[1])

            if id_1 >= total_nodes:
                nodes_set_1 = clusters_collection[id_1 - total_nodes]
            else:
                nodes_set_1 = [idx_2_id[id_1]]

            if id_2 >= total_nodes:
                nodes_set_2 = clusters_collection[id_2 - total_nodes]
            else:
                nodes_set_2 = [idx_2_id[id_2]]

            clusters_collection.append(nodes_set_1 + nodes_set_2)

        # what each cluster contains, average information flow between nodes and nodes in cluster
        return clusters_collection, 1./clust_linkmap[:, 2], clust_linkmap[:, 3]