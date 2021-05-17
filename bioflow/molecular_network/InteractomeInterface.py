"""
This module contains all the routines that are responsible for pulling
the matrices out of the knowledge time and routines on them
"""
# from pympler import muppy, summary, tracker

import hashlib
import itertools
import json
import os
import pickle
import string
from collections import defaultdict
from copy import copy
from random import shuffle, sample
from time import time
import numpy as np
from scipy.sparse import lil_matrix
import scipy.sparse as spmat
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from itertools import chain
from typing import Union, Tuple, List

import bioflow.configs.main_configs as confs
from bioflow.utils.gdfExportInterface import GdfExportInterface
from bioflow.utils.io_routines import write_to_csv, dump_object, undump_object
from bioflow.utils.log_behavior import get_logger

from bioflow.sample_storage.mongodb import insert_interactome_rand_samp
from bioflow.configs.main_configs import internal_storage
from bioflow.algorithms_bank import conduction_routines as cr
from bioflow.algorithms_bank import weigting_policies as wp
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.algorithms_bank import sampling_policies
from bioflow.algorithms_bank.flow_calculation_methods import general_flow,\
    reduce_and_deduplicate_sample, evaluate_ops, reduce_ops
from bioflow.algorithms_bank.sampling_policies import characterize_flow_parameters, _is_int


log = get_logger(__name__)


class InteractomeInterface(object):
    """
    Interface between interactome in the knowledge database and the interactome graph laplacian

    :param background_up_ids: (optional) background that will be used for sparse_sampling of random
    nodes to build a comparison interface for the
    """

    def __init__(self, background_up_ids=()): # TRACING: [weighted background] GOOD
        self.init_time = time()
        self.partial_time = time()

        self.adjacency_matrix = np.zeros((4, 4))
        # This is just non-normalized laplacian matrix
        self.laplacian_matrix = np.zeros((4, 4))
        self.non_norm_laplacian_matrix = np.zeros((4, 4))

        # REFACTOR [stateless]: decouple into argument
        self.adj_eigenvects = np.zeros((4, 4))
        self.adj_eigenvals = np.zeros((4, 4))
        self.cond_eigenvects = np.zeros((4, 4))
        self.cond_eigenvals = np.zeros((4, 4))

        self.neo4j_id_2_matrix_index = {}
        self.matrix_index_2_neo4j_id = {}
        self.neo4j_id_2_display_name = {}
        self.neo4j_id_2_legacy_id = {}
        self.neo4j_id_2_node_type = {}
        self.neo4j_id_2_localization = {}
        self.known_uniprots_neo4j_ids = []

        self.maps_dumps_location = confs.Dumps.interactome_maps
        self.adjacency_dumps_location = confs.Dumps.interactome_adjacency_matrix
        self.laplacian_dumps_location = confs.Dumps.interactome_laplacian_matrix

        char_set = string.ascii_uppercase + string.digits
        self.thread_hex = ''.join(sample(char_set * 6, 6))

        self.UP2UP_voltages = {}
        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}

        # REFACTOR [stateless]: decouple into arguments
        self._active_up_sample: List[int] = []
        self._active_weighted_sample: List[Tuple[int, float]] = []
        self._secondary_weighted_sample: Union[None, List[Tuple[int, float]]] = None

        self._flow_calculation_method = general_flow
        self._ops_evaluation_method = evaluate_ops
        self._ops_reduction_method = reduce_ops

        self.sparsely_sampled = False  # TRACING: sparse sampling save location.
        # TRACING: always false, would have been used for export otherwise
        self._background = background_up_ids  # TRACING: [weighted background] GOOD
        log.debug('_background set to %d' % len(background_up_ids))

    def pretty_time(self):
        """
        Times the execution

        :return: tuple containing the time since the creation of the Matrix_getter object and
         since the last cal of function formatted as string
        """
        init_time, partial_time = (round(time() - self.init_time),
                                   round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (
            int(init_time) // 60, init_time % 60, int(partial_time) // 60, partial_time % 60)
        self.partial_time = time()
        return pload

    def _time(self):
        """
        time since previous partial timer flag was reset.
        """
        partial_time = time() - self.partial_time
        return partial_time

    def reset_thread_hex(self):
        """
        Reset the hex identifier of the object. used when multiprocessing

        :return:
        """
        char_set = string.ascii_uppercase + string.digits
        self.thread_hex = ''.join(sample(char_set * 6, 6))

    # REFACTOR: [better environment]:  potential improvement for all dump/undump methods
    #  is to use mongo storage and index from environment not to re-compute them every time
    def _dump_matrices(self):
        """
        dumps self.adjacency_matrix and self.laplacian_matrix
        """
        dump_object(confs.Dumps.interactome_adjacency_matrix, self.adjacency_matrix)
        dump_object(confs.Dumps.interactome_laplacian_matrix, self.laplacian_matrix)

    def _undump_matrices(self):
        """
        undumps self.adjacency_matrix and self.laplacian_matrix
        """
        self.adjacency_matrix = undump_object(confs.Dumps.interactome_adjacency_matrix)
        self.laplacian_matrix = undump_object(confs.Dumps.interactome_laplacian_matrix)
        self.non_norm_laplacian_matrix = self.laplacian_matrix.copy()

    def _dump_eigen(self):
        """
        dumps self.adj_eigenvals and self.laplacian_matrix and writes them to csv
        """
        write_to_csv(confs.Dumps.eigen_VaMat, self.adj_eigenvals)
        write_to_csv(confs.Dumps.eigen_ConMat, self.cond_eigenvals)
        dump_object(confs.Dumps.val_eigen, (self.adj_eigenvals, self.adj_eigenvects))
        dump_object(
            confs.Dumps.cond_eigen,
            (self.cond_eigenvals,
             self.cond_eigenvects))

    def _undump_eigen(self):
        """
        undumps self.adj_eigenvals and self.laplacian_matrix
        """
        self.adj_eigenvals, self.adj_eigenvects = undump_object(
            confs.Dumps.val_eigen)
        self.cond_eigenvals, self.cond_eigenvects = undump_object(
            confs.Dumps.cond_eigen)

    def _dump_maps(self):
        """
        dumps all the elements required for the mapping between the types and ids
         of database entries and matrix columns
        """
        log.debug("pre-dump e_p_u_b_i length: %s", len(self._active_up_sample))
        log.debug("dumping into: %s", confs.Dumps.interactome_maps)
        dump_object(
            confs.Dumps.interactome_maps,
            (self.neo4j_id_2_matrix_index,
             self.matrix_index_2_neo4j_id,
             self.neo4j_id_2_display_name,
             self.neo4j_id_2_legacy_id,
             self.neo4j_id_2_node_type,
             self.neo4j_id_2_localization,
             self.known_uniprots_neo4j_ids,
             self._active_up_sample))

    def _undump_maps(self):
        """
        undumps all the elements required for the mapping between the types and ids of
        database entries and matrix columns
        """
        log.debug("undumping from %s", confs.Dumps.interactome_maps)
        self.neo4j_id_2_matrix_index, self.matrix_index_2_neo4j_id, \
        self.neo4j_id_2_display_name, self.neo4j_id_2_legacy_id, self.neo4j_id_2_node_type, \
        self.neo4j_id_2_localization, \
        self.known_uniprots_neo4j_ids, self._active_up_sample = \
            undump_object(confs.Dumps.interactome_maps)
        log.debug("post-undump e_p_u_b_i length: %s", len(self._active_up_sample))

    def _dump_memoized(self):    # TRACING: add weighted samples chars
        """
        In a JSON, stores a dump of the following properties:
        {
            'UP_hash': md5,
            'sys_hash': self.md5_hash(),
            'size': len(self._active_up_sample),
            'UPs': pickle.dumps(self._active_up_sample),
            'currents': pickle.dumps((self.current_accumulator, self.node_current))}
        :return:
        """
        md5 = hashlib.md5(
            json.dumps(
                sorted(
                    self._active_up_sample),
                sort_keys=True).encode('utf-8')).hexdigest()

        payload = {
            'UP_hash': md5,
            'sys_hash': self.md5_hash(),
            'size': len(self._active_up_sample),
            'UPs': pickle.dumps(self._active_up_sample),
            'currents': pickle.dumps((self.current_accumulator, self.node_current))}
        dump_object(confs.Dumps.Interactome_Analysis_memoized, payload)

    @staticmethod
    def _undump_memoized() -> dict:    # TRACING: add weighted samples chars
        """
        Retrieves a JSON dump of the following properties:
        {
            'UP_hash': md5,
            'sys_hash': self.md5_hash(),
            'size': len(self._active_up_sample),
            'UPs': pickle.dumps(self._active_up_sample),
            'currents': pickle.dumps((self.current_accumulator, self.node_current)),
            'voltages': pickle.dumps(self.uniprots_2_voltage)}

        :return:
        """
        return undump_object(confs.Dumps.Interactome_Analysis_memoized)

    def normalize_laplacian(self):
        """
        Performs Laplacian Normalization

        :return:
        """
        self.non_norm_laplacian_matrix = self.laplacian_matrix.copy()
        D = self.laplacian_matrix.diagonal()
        mD = np.power(D, -0.5)
        lpl_shape = self.laplacian_matrix.shape
        D = lil_matrix(lpl_shape)
        D.setdiag(mD)
        D = D.tocsc()
        self.laplacian_matrix = (D.dot(self.laplacian_matrix)).dot(D)

    def create_val_matrix(self,
                          node_dict: dict, edge_list: list,
                          adj_weight_policy_function=wp.active_default_adj_weighting_policy,
                          lapl_weight_policy_function=wp.active_default_lapl_weighting_policy):
        """
        Creates the ajacency/laplacian matrix pair for the given set of nodes and edges and the
        weighting policy functions

        :param node_dict: dictionary mapping node id to node objects
        :param edge_list: list of edges that exist in the database
        :param adj_weight_policy_function: adjacency matrix weight policy function
        :param lapl_weight_policy_function: laplacian matrix weight policy function
        :return: node_id to matrix idx map, matrix_idx to node_id map, adjacency matrix,
        laplacian matrix
        """
        log.debug('new val matrix: %dx%d' % (len(node_dict), len(node_dict)))
        adjacency_matrix = lil_matrix((len(node_dict), len(node_dict)), dtype=np.float)
        laplacian_matrix = lil_matrix((len(node_dict), len(node_dict)), dtype=np.float)

        node_id_2_mat_idx = {_id: _i for _i, _id in enumerate(node_dict.keys())}
        mat_idx_2_note_id = {_i: _id for _id, _i in node_id_2_mat_idx.items()}

        for rel_obj in edge_list:
            # we need to use ids because the node objects stored by the node are property-free
            from_id = rel_obj.start_node.id
            to_id = rel_obj.end_node.id

            from_node = node_dict[from_id]
            to_node = node_dict[to_id]

            lapl_weight = lapl_weight_policy_function(from_node, to_node, rel_obj)
            adj_weight = adj_weight_policy_function(from_node, to_node, rel_obj)

            from_idx = node_id_2_mat_idx[from_id]
            to_idx = node_id_2_mat_idx[to_id]

            adjacency_matrix[from_idx, to_idx] += adj_weight
            adjacency_matrix[to_idx, from_idx] += adj_weight

            laplacian_matrix[from_idx, to_idx] -= lapl_weight
            laplacian_matrix[to_idx, from_idx] -= lapl_weight
            laplacian_matrix[from_idx, from_idx] += lapl_weight
            laplacian_matrix[to_idx, to_idx] += lapl_weight

        return node_id_2_mat_idx, mat_idx_2_note_id, adjacency_matrix, laplacian_matrix

    def giant_component_node_idxs(self, adjacency_matrix):
        """
        Finds the indexes of the lines and columns in the adjacency matrix that are in the giant
        component

        :param adjacency_matrix: the adjacency matrix where to find the giant component
        :return:
        """
        component_ids, node_2_id = connected_components(adjacency_matrix, directed=False)

        counters = np.zeros((component_ids,))

        for id in range(0, component_ids):
            counters[id] = np.sum((node_2_id == id).astype(int))

        biggest_component_id = np.argmax(counters)

        idx_in_g_component = np.arange(adjacency_matrix.shape[0])[node_2_id == biggest_component_id]
        
        return idx_in_g_component.tolist()

    def full_rebuild(self, adj_weight_policy_function=wp.active_default_adj_weighting_policy,
                           lapl_weight_policy_function=wp.active_default_lapl_weighting_policy):
        """
        Performs a complete rebuild of the InterfaceClass Instance based on parameters provided
        upon construction based on the data in the knowledge database. Upon rebuild saves a copy
        that can be rapidly resurrected with the fast_load() method

        :param adj_weight_policy_function: adjacency matrix weight policy function
        :param lapl_weight_policy_function: laplacian matrix weight policy function

        :return: None
        """
        # giant component recomputation and writing
        DatabaseGraph.erase_node_properties(['main_connex'])

        all_nodes_dict, edges_list = DatabaseGraph.parse_physical_entity_net(main_connex_only=False)

        _, mat_idx_2_note_id, adjacency_matrix, _ = \
            self.create_val_matrix(all_nodes_dict, edges_list, wp.flat_policy, wp.flat_policy)

        giant_component_mat_indexes = self.giant_component_node_idxs(adjacency_matrix)

        giant_component_db_ids = [mat_idx_2_note_id[_idx] for _idx in giant_component_mat_indexes]

        DatabaseGraph.batch_set_attributes(giant_component_db_ids,
                                           [{'main_connex': 'True'}]*len(giant_component_db_ids))

        # only giant component parsing
        nodes_dict, edges_list = DatabaseGraph.parse_physical_entity_net(main_connex_only=True)

        node_id_2_mat_idx, mat_idx_2_note_id, adjacency_matrix, laplacian_matrix = \
            self.create_val_matrix(nodes_dict, edges_list,
                                   adj_weight_policy_function=adj_weight_policy_function,
                                   lapl_weight_policy_function=lapl_weight_policy_function)

        all_uniprot_nodes = DatabaseGraph.get_all('UNIPROT')

        self.adjacency_matrix = adjacency_matrix
        self.laplacian_matrix = laplacian_matrix

        self.get_eigen_spectrum(100)

        self.neo4j_id_2_matrix_index = node_id_2_mat_idx
        self.matrix_index_2_neo4j_id = mat_idx_2_note_id

        self.neo4j_id_2_display_name = {_id: _node['displayName']
                                        for (_id, _node) in nodes_dict.items()}
        self.neo4j_id_2_legacy_id = {_id: _node['legacyID']
                                     for (_id, _node) in nodes_dict.items()}
        self.neo4j_id_2_node_type = {_id: list(_node.labels)[0]
                                     for (_id, _node) in nodes_dict.items()}
        self.neo4j_id_2_localization = {_id: _node.get('localization', 'NA')
                                        for (_id, _node) in nodes_dict.items()}

        self.known_uniprots_neo4j_ids = [_node_id for _node_id
                                         in self.neo4j_id_2_matrix_index.keys()]

        self._active_up_sample = []  # REFACTOR [stateless]: decouple into argument

        self._dump_maps()  # DONE
        self._dump_matrices()  # DONE
        self._dump_eigen()  # DONE

    def get_eigen_spectrum(self, biggest_eigvals_to_get):
        """
        Recovers the eigenspectrum associated to the *n* biggest eigenvalues, where *n* is
        specified by biggest_eigvals_to_get. If the Adjacency and conductance matrix haven't
        been preloaded first, will raise an Exception

        :param biggest_eigvals_to_get: specifies how many biggest eigenvalues we are willing to get.
        :raise Exception: "Matrix must be pre-loaded first" if self.adjacency_matrix and
        self.laplacian_matrix have not been computed anew or pre-loaded first
        """
        if self.laplacian_matrix.shape == (4, 4):
            log.critical("Matrix must be pre-loaded first")
            raise Exception("Matrix must be pre-loaded first")

        log.info("entering eigenvect computation; %s", self.pretty_time())

        self.adj_eigenvals, self.adj_eigenvects = eigsh(
            self.adjacency_matrix, biggest_eigvals_to_get)
        self.cond_eigenvals, self.cond_eigenvects = eigsh(
            self.laplacian_matrix, biggest_eigvals_to_get)

        log.debug("Adjacency matrix eigenvalues:")
        log.debug(self.adj_eigenvals)
        log.debug('<======================>')
        log.debug("Laplacian matrix eigenvalues:")
        log.debug(self.cond_eigenvals)
        log.debug('<======================>')
        log.debug("all laplacian eigenvalues above%s", np.all(eigsh(self.laplacian_matrix)[0] > 0))
        log.info("Finished eigenvalues computation, starting the dump %s", self.pretty_time())

    def fast_load(self):
        """
        Rapidly resurrects the InterfaceClass Instance based on parameters provided
        upon construction. If parameters are mismatched, raises exceptions signalling what
        parameters were mismatched. Trims the background provided upon construction down to what
        actually be sampled (self._background)

        """
        self._undump_maps()
        self._undump_matrices()
        self._undump_eigen()

        if self._background:  # TRACING: [weighted background] GOOD
            if _is_int(self._background[0]):
                self._background = list(set(self.known_uniprots_neo4j_ids).intersection(
                    set(self._background)))
            else:
                self._background = [(_id, _weight)
                                    for _id, _weight in self._background
                                    if _id in self.known_uniprots_neo4j_ids]

        else:  # TRACING: [weighted background] GOOD
            self._background = list(set(self.known_uniprots_neo4j_ids))

    def get_descriptor_for_index(self, index):
        """
        Unused.

        Recovers a descriptor set for a given index in the current matrix mapping

        :param index: idenx for which return a desciptor
        :return: Type, displayName and if a localization is given, returns display name too.
        :rtype: tuple
        """
        if self.matrix_index_2_neo4j_id[index] in list(self.neo4j_id_2_localization.keys()):
            return (self.neo4j_id_2_node_type[self.matrix_index_2_neo4j_id[index]],
                    self.neo4j_id_2_display_name[self.matrix_index_2_neo4j_id[index]],
                    self.neo4j_id_2_localization[self.matrix_index_2_neo4j_id[index]])
        else:
            return (self.neo4j_id_2_node_type[self.matrix_index_2_neo4j_id[index]],
                    self.neo4j_id_2_display_name[self.matrix_index_2_neo4j_id[index]])

    def md5_hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been correctly
        defined before dump/retrieval
        """
        sorted_initial_set = sorted(self.neo4j_id_2_matrix_index.keys())

        # REFACTOR: should involve database used for the build metadata
        data = [
            sorted_initial_set,
            confs.line_loss,
            confs.use_normalized_laplacian,
            confs.fraction_edges_dropped_in_laplacian,
            (confs.env_skip_reactome, confs.env_skip_hint, confs.env_skip_biogrid)]

        md5 = hashlib.md5(json.dumps(data, sort_keys=True).encode('utf-8')).hexdigest()

        log.debug("System md5 hashing done. hash: %s. parameters: \n"
                  "\t sorted init set hash: %s \n"
                  "\t line loss: %s\n"
                  "\t l_norm: %s\n"
                  "\t edge drop: %s\n"
                  "\t skipping reactome|hint|biogrid: %s\n"
                  % (md5,
                     hashlib.md5(json.dumps(sorted_initial_set, sort_keys=True).encode('utf-8')).hexdigest(),
                     confs.line_loss,
                     confs.use_normalized_laplacian,
                     confs.fraction_edges_dropped_in_laplacian,
                     (confs.env_skip_reactome, confs.env_skip_hint, confs.env_skip_biogrid)))

        return str(md5)

    def active_sample_md5_hash(self, sparse_rounds):
        sys_hash = self.md5_hash()

        if self._background: # TRACING: [weighted background]: GOOD
            if _is_int(self._background[0]):
                background = sorted(self._background)
            else:
                background = sorted(self._background, key=lambda x: x[1])

        sample_chars = characterize_flow_parameters(self._active_weighted_sample,
                                                    self._secondary_weighted_sample,
                                                    sparse_rounds)

        hashlib.md5(json.dumps(background, sort_keys=True).encode(
                                                'utf-8')).hexdigest()

        data = [
            sys_hash,
            background,
            self._flow_calculation_method.__name__,
            sparse_rounds,
            sample_chars[0], sample_chars[1], sample_chars[2],
            sample_chars[3], sample_chars[4], sample_chars[5],
            sample_chars[6], sample_chars[7]
        ]

        md5 = hashlib.md5(json.dumps(data, sort_keys=True).encode('utf-8')).hexdigest()

        log.debug('Active sample md5 hashing done: %s. parameters: \n'
                  '\tsys_hash: %s\n'
                  '\tbackground hash: %s\n'
                  '\tflow policy: %s\n'
                  '\tsparse sampling: %s\n'
                  '\tmain sample chars: %d/%d/%s\n'
                  '\tsec sample chars: %d/%d/%s\n'
                  '\tsparse sampling: %d\n'
                  '\toverall hash: %s' % (md5,
                                          sys_hash,
                                          hashlib.md5(
                                              json.dumps(background, sort_keys=True).encode(
                                                'utf-8')).hexdigest(),
                                          self._flow_calculation_method.__name__,
                                          sparse_rounds,
                                          sample_chars[0], sample_chars[1], sample_chars[2],
                                          sample_chars[3], sample_chars[4], sample_chars[5],
                                          sample_chars[6], sample_chars[7]))

        return str(md5)

    def set_flow_sources(self, sample, secondary_sample):

        def _verify_uniprot_ids(uniprot_vector: List[Tuple[int, float]]):
            # TRACING: rename to id_weight_vector
            uniprots = np.array(uniprot_vector)[:, 0].tolist()

            if not set(uniprots) <= set(self.known_uniprots_neo4j_ids):

                log.warn('Following reached uniprots neo4j_ids were not retrieved upon the '
                         'circulation matrix construction: \n %s',
                         (set(uniprots) - set(self.known_uniprots_neo4j_ids)))

            _filter = [True
                       if uniprot in self.known_uniprots_neo4j_ids
                       else False
                       for uniprot in uniprots]

            return np.array(uniprot_vector)[_filter, :].tolist()

        self._active_weighted_sample = _verify_uniprot_ids(reduce_and_deduplicate_sample(sample))

        self._active_up_sample = np.array(self._active_weighted_sample)[:, 0].tolist()

        if secondary_sample is not None:
            self._secondary_weighted_sample = \
                _verify_uniprot_ids(reduce_and_deduplicate_sample(secondary_sample))

            self._active_up_sample = list(set(self._active_up_sample
                                              + np.array(self._secondary_weighted_sample)[:, 0].tolist()))

    def evaluate_ops(self, sparse_rounds=-1):
        log.debug('evaluate_ops call')
        ro = sampling_policies.characterize_flow_parameters(self._active_weighted_sample,
                                                            self._secondary_weighted_sample,
                                                            -2)
        return self._ops_evaluation_method(ro[1], ro[3], sparse_rounds)

    def reduce_ops(self, ops_limit):
        log.debug('reduce_ops call')
        ro = sampling_policies.characterize_flow_parameters(self._active_weighted_sample,
                                                            self._secondary_weighted_sample,
                                                            -2)

        return self._ops_reduction_method(ro[1], ro[3], ops_limit)

    def set_uniprot_source(self, uniprots):
        """
        Sets the deprecated_reached_uniprots_neo4j_id_list on which the circulation computation routines will
        be performed by the other methods. Avoids passing as argument large lists of parameters.

        :param uniprots: List of node IDs of the uniprots on which we would like to
        perform current computations
        :raise Warning: if the uniprots were not present in the set of GOs for which
        we built the system or had no GO attached to them
        """
        if not set(uniprots) <= set(self.known_uniprots_neo4j_ids):

            log.warn('Following reached uniprots neo4j_ids were not retrieved upon the '
                     'circulation matrix construction: \n %s',
                     (set(uniprots) - set(self.known_uniprots_neo4j_ids)))

        self._active_up_sample = \
            [uniprot for uniprot in uniprots if uniprot in self.known_uniprots_neo4j_ids]

    def compute_current_and_potentials(
            self,
            memoized: bool = True,  # is required to enable the fast loading.
            incremental: bool = False,  # This is always false and was used in order to resume the
            # sparse_sampling
            cancellation: bool = True,
            sparse_samples: int = -1,
            fast_load: bool = False):  # REFACTOR: this should not be implemented
        # this way.
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier
        knowledge flow analysis

        :param memoized: if the tensions between individual nodes and voltages will be
            remembered - required for clustering. Incompatible with `sparse_sample=True`
        :param incremental: if True, all the circulation computation will be added to the
            existing ones. Useful for the computation of particularly big systems with
            intermediate dumps
        :param cancellation: divides the final current by number of bioflow-sink pairs
        :param sparse_samples: if set to an integer >0, the sparse_sampling will be sparse and not
            dense,i.e. instead of computation for each node pair, only an estimation will be made,
            equal to computing sparse_samples association with other randomly chosen nodes
        :param fast_load: if True, will try to lad a pre-saved instance
        :return: adjusted conduction system
        """

        if confs.use_normalized_laplacian:
            self.normalize_laplacian()

        if fast_load:
            payload = self._undump_memoized()  # TRACING: add weighted samples chars

              # TRACING: hash of samples. Fast ressurection requires self._weighted_samples
            UP_hash = hashlib.md5(
                json.dumps(
                    sorted(self._active_up_sample),
                    sort_keys=True).encode("utf-8")
                ).hexdigest()

            if payload['sys_hash'] == self.md5_hash() and payload['UP_hash'] == UP_hash:
                self.current_accumulator, self.node_current = pickle.loads(payload['currents'])

            index_current = cr.get_current_through_nodes(self.current_accumulator)
            log.info('current accumulator shape %s', self.current_accumulator.shape)

            self.node_current.update(
                dict((self.matrix_index_2_neo4j_id[idx], val) for idx, val in
                     enumerate(index_current)))

            return None

        if not incremental or self.current_accumulator == np.zeros((2, 2)):
            self.current_accumulator = spmat.csc_matrix(self.laplacian_matrix.shape)
            self.UP2UP_voltages = {}
            self.node_current = defaultdict(float)

        translated_active_weighted_sample = [(self.neo4j_id_2_matrix_index[UP], _w) for
                                             UP, _w in self._active_weighted_sample]

        if self._secondary_weighted_sample is not None:
            translated_secondary_weighted_sample = [(self.neo4j_id_2_matrix_index[UP], _w) for
                                                    UP, _w in self._secondary_weighted_sample]

        else:
            translated_secondary_weighted_sample = None

        current_accumulator, up_pair_2_voltage = \
            cr.main_flow_calc_loop(self.laplacian_matrix,
                                   translated_active_weighted_sample,
                                   secondary_sample=translated_secondary_weighted_sample,
                                   cancellation=cancellation,
                                   sparse_rounds=sparse_samples,
                                   potential_diffs_remembered=True,
                                   thread_hex=self.thread_hex,
                                   active_sampling_function=self._flow_calculation_method)

        self.UP2UP_voltages.update(
                dict(((self.matrix_index_2_neo4j_id[i],
                       self.matrix_index_2_neo4j_id[j]),
                      voltage)
                     for (i, j), voltage in up_pair_2_voltage.items()))

        # if sparse_samples > 1:
        #     current_accumulator, _ = \
        #         cr.main_flow_calc_loop(self.laplacian_matrix,
        #                                [self.neo4j_id_2_matrix_index[UP]
        #                                 for UP in self._active_up_sample],
        #                                cancellation=cancellation,
        #                                sparse_rounds=sparse_samples,
        #                                thread_hex=self.thread_hex)
        #
        # else:
        #     current_accumulator, up_pair_2_voltage =\
        #         cr.main_flow_calc_loop(self.laplacian_matrix,
        #                                [self.neo4j_id_2_matrix_index[UP]
        #                                 for UP in self._active_up_sample],
        #                                cancellation=cancellation,
        #                                potential_diffs_remembered=True,
        #                                thread_hex=self.thread_hex)
        #
        #
        #     self.UP2UP_voltages.update(
        #         dict(((self.matrix_index_2_neo4j_id[i],
        #                self.matrix_index_2_neo4j_id[j]),
        #               voltage)
        #              for (i, j), voltage in up_pair_2_voltage.items()))

        if incremental:
            self.current_accumulator = self.current_accumulator + current_accumulator
        else:
            self.current_accumulator = current_accumulator

        index_current = cr.get_current_through_nodes(self.current_accumulator)
        log.info('current accumulator shape %s, sum %s', current_accumulator.shape, np.sum(current_accumulator))
        self.node_current.update(
            dict((self.matrix_index_2_neo4j_id[idx], val) for idx, val in enumerate(index_current)))

        if memoized:
            self._dump_memoized()

    def format_node_props(self, node_current, limit=0.01):
        """
        Formats the nodes for the analysis by in the knowledge_access_analysis module

        :param node_current: Current through the entity
        :param limit: hard limit to go_namespace_filter out the GO terms with too little current
        (compensates the minor currents in the gird)
        :return: {Entity:[node current, node degree]}
        """
        characterization_dict = {}
        limit_current = max(node_current.values()) * limit
        for NodeID, i in self.neo4j_id_2_matrix_index.items():
            if node_current[NodeID] > limit_current:
                characterization_dict[NodeID] = [node_current[NodeID],
                                                 self.non_norm_laplacian_matrix[i, i]]
        return characterization_dict

    def export_conduction_system(self,
                                 p_value_dict: dict = None,
                                 output_location: str = ''):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format and
         flushes it into a file that can be viewed with Gephi

         :param p_value_dict:
         :param output_location:
        """

        if self.sparsely_sampled:
            log.warning('Computation of the information circulation was not complete, %s',
                        'most likely due to the sparse_sampling')

        node_char_names = [
            'Current',
            'Type',
            'Legacy_ID',
            'Names',
            'Degree',
            'Source',
            'p-value',
            'p_p-value',
            'rel_value',
            'std_diffs']

        node_char_types = [
            'DOUBLE',
            'VARCHAR',
            'VARCHAR',
            'VARCHAR',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE']  # TRACING: add the weights at the start

        if p_value_dict is None:
            p_value_dict = defaultdict(lambda: (np.nan, np.nan, np.nan))

        nan_neg_log10 = lambda x: x if type(x) is str else -np.log10(x)

        characterization_dict = {}

        log.info("Laplacian shape %s", self.laplacian_matrix.shape)
        log.info("Matrix size %s", max(self.matrix_index_2_neo4j_id.keys()))

        for NodeID in self.node_current.keys():
            matrix_index = self.neo4j_id_2_matrix_index[NodeID]

            if NodeID not in list(self.neo4j_id_2_display_name.keys()):
                log.warning('neo4j id %s does not seem to appear in the main import set', NodeID)
                log.warning('corresponding matrix id is %s', matrix_index)
                continue

            if not self.neo4j_id_2_display_name[NodeID]:
                log.warning('neo4j id %s maps to a display ID that is void', NodeID)
                log.warning('corresponding matrix id is %s', matrix_index)
                self.neo4j_id_2_display_name[NodeID] = "None"
                continue

            characterization_dict[NodeID] = [
                str(self.node_current[NodeID]),
                self.neo4j_id_2_node_type[NodeID],
                self.neo4j_id_2_legacy_id[NodeID],
                self.neo4j_id_2_display_name[NodeID].replace(',', '-'),
                str(self.laplacian_matrix[matrix_index, matrix_index]),
                str(float(int(NodeID in self._active_up_sample))), # TRACING: add the weights at the start
                # TRACING: add if it is a secondary sample or not.
                # CURRENTPASS: just map positive from primary_sample and negative from sec_sample
                str(p_value_dict[int(NodeID)][0]),
                str(nan_neg_log10(p_value_dict[int(NodeID)][0])),
                str(float(p_value_dict[int(NodeID)][1])),
                str(float(p_value_dict[int(NodeID)][2]))]

        if output_location == '':
            output_location = confs.NewOutputs().Interactome_GDF_output


        gdf_exporter = GdfExportInterface(
            target_fname=output_location,
            field_names=node_char_names,
            field_types=node_char_types,
            node_properties_dict=characterization_dict,
            min_current=0.0001,
            index_2_label=self.matrix_index_2_neo4j_id,
            label_2_index=self.neo4j_id_2_matrix_index,
            current_matrix=self.current_accumulator)
            # TODO: [Better stats]: twister compared to random sample?
        gdf_exporter.write()

    def randomly_sample(
            self,
            iterations,
            sparse_rounds=-1,
            no_add=False,
            pool_no=None,
            sampling_policy=sampling_policies.matched_sampling,
            optional_sampling_param = 'exact'):
        """
        Randomly samples the set of deprecated_reached_uniprots_neo4j_id_list used to create the model.
        This is the null model creation routine. It will match the null model to the currently
        loaded active samples.

        :param iterations: how many times we would like to sample each uniprot number
        :param sparse_rounds:  if we want to use sparse sparse_sampling (useful in case of
        large uniprot sets), we would use this option
        :param no_add: if set to True, the result of sparse_sampling will not be added to the
        database of samples. Useful if re-running tests with similar parameters several times.
        :param pool_no: explicit sparse_sampling pool number (used for reporting/debugging)
        :param sampling_policy: sampling policy used
        :param sampling_policy_options: sampling policy optional argument
        :raise Exception: if the number of items in the samples size ann samples_each size
        are different
        """
        sample_size = -1 # for legacy reasons

        # TODO: [Better sparse_sampling]: include the limitations on the types of nodes to sample:

        sample_chars = characterize_flow_parameters(self._active_weighted_sample,
                                                    self._secondary_weighted_sample, sparse_rounds)

        super_hash = sample_chars[7]

        log.info('Starting a random sampler: \n'
                 '\tsampling policy & optional param: %s/%s\n'
                 '\tflow policy: %s\n'
                 '\tsparse sampling: %s\n'
                 '\tmain sample chars: %d/%d/%s\n'
                 '\tsec sample chars: %d/%d/%s\n'
                 '\tsparse sampling: %d\t'
                 '\toverall hash: %s' % (sampling_policy.__name__, optional_sampling_param,
                                         self._flow_calculation_method.__name__,
                                         sparse_rounds,
                                         sample_chars[0], sample_chars[1], sample_chars[2],
                                         sample_chars[3], sample_chars[4], sample_chars[5],
                                         sample_chars[6], sample_chars[7]))

        preserved_sample = self._active_weighted_sample.copy()

        if self._secondary_weighted_sample is not None:
            preserved_sec_sample = self._secondary_weighted_sample.copy()

        else:
            preserved_sec_sample = None

        for i, sample, sec_sample in sampling_policy(preserved_sample,
                                                     preserved_sec_sample,
                                                     # TRACING: [weighted background] GOOD
                                                     self._background,
                                                     iterations,
                                                     optional_sampling_param):

            # print('debug: selected UProt IDs :', sample)

            self.set_flow_sources(sample, sec_sample)

            sample_chars = characterize_flow_parameters(sample, sec_sample, sparse_rounds)
            sample_hash = sample_chars[-1]

            # self.set_uniprot_source(sample)

            log.info('Sampling thread: %s, Thread hex: %s; Random sample %d/%d \n'
                     'sparse_sampling characteristics: sys_hash: %s, sample_hash: %s, '
                     'target_hash: %s' %
                     (pool_no, self.thread_hex, i, iterations,
                      self.md5_hash(), sample_hash, super_hash))

            # TODO: [load bar]: the external loop progress bar goes here

            self.compute_current_and_potentials(memoized=False, sparse_samples=sparse_rounds)

            sample_ids_md5 = hashlib.md5(
                json.dumps(
                    sorted(self._active_up_sample),
                    sort_keys=True).encode('utf-8')).hexdigest()

            if not no_add:
                log.info("Sampling thread %s: Adding a blanc:"
                         "\t sys hash: %s "
                         "\t sample hash: %s \t active sample hash: %s \t target_hash: %s \t "
                         "sparse_rounds: %s \t sampling policy: %s\t sampling_options: %s \t "
                         "matrix weight: %s"
                         % (pool_no, self.md5_hash(),
                            sample_hash, self.active_sample_md5_hash(sparse_rounds), super_hash,
                            sparse_rounds, sampling_policy.__name__, optional_sampling_param,
                            np.sum(self.current_accumulator)))

                insert_interactome_rand_samp(  # INTEST: sample storage change
                    {
                        'UP_hash': sample_ids_md5,  # specific retrieval, but inexact.
                        'sys_hash': self.md5_hash(),
                        'active_sample_hash': self.active_sample_md5_hash(sparse_rounds),
                        'target_sample_hash': super_hash,
                        'sampling_policy': sampling_policy.__name__,
                        'sampling_policy_options': optional_sampling_param,
                        'size': -1,  # TRACING: to be removed
                        'sparse_rounds': sparse_rounds,
                        'UPs': pickle.dumps(self._active_up_sample),
                        'sample': pickle.dumps(self._active_weighted_sample),
                        'sec_sample': pickle.dumps(self._secondary_weighted_sample),
                        'currents': pickle.dumps(
                            (self.current_accumulator,
                             self.node_current)),  # TRACING: node currents are dead: deprecate
                        'voltages': pickle.dumps(
                            self.UP2UP_voltages)})

            # if not sparse_rounds:
            #     log.info('Sampling thread %s: Thread hex: %s \t Sample size: %s \t iteration: %s\t compop/s: %s \t '
            #              'time: %s ',
            #              pool_no, self.thread_hex, sample_size, i,
            #              "{0:.2f}".format(sample_size * (sample_size - 1) / 2 / self._time()),
            #              self.pretty_time())
            # else:
            #     log.info('Sampling thread %s: Thread hex: %s \t Sample size: %s \t iteration: %.2f \t compop/s: %s \t '
            #              'time: %s, sparse @ %s ',
            #              pool_no, self.thread_hex, sample_size, i,
            #              "{0:.2f}".format(sample_size * sparse_rounds / 2 / self._time()),
            #              self.pretty_time(), sparse_rounds)

        self._active_weighted_sample = preserved_sample
        self._secondary_weighted_sample = preserved_sec_sample



    @staticmethod
    def compare_dumps(dumps_folder_1, dumps_folder_2):
        """
        Compares the dumps between two folders (debug function for methods upgrade)

        :param dumps_folder_1: old dumps
        :param dumps_folder_2: new dumps
        :return:
        """
        neo4j_id_2_matrix_index_1,\
        matrix_index_2_neo4j_id_1, _,\
        neo4j_id_2_legacy_id_1, _, _,\
        reached_uniprots_neo4j_id_list_1,\
        all_uniprots_neo4j_id_list_1, _, _, _, _,\
        _ = undump_object(dumps_folder_1+'/dump2.dump')

        adjacency_matrix_1 = undump_object(dumps_folder_1+'/pickleDump3.dump')
        legacy_id_2_neo4j_id_1 = dict((value, key) for key, value in neo4j_id_2_legacy_id_1.items())

        neo4j_id_2_matrix_index_2, \
        matrix_index_2_neo4j_id_2, _, \
        neo4j_id_2_legacy_id_2, _, _, \
        reached_uniprots_neo4j_id_list_2, \
        all_uniprots_neo4j_id_list_2, _, _, _, _, \
        _ = undump_object(dumps_folder_2 + '/dump2.dump')

        adjacency_matrix_2 = undump_object(dumps_folder_2 + '/pickleDump3.dump')
        legacy_id_2_neo4j_id_2 = dict((value, key) for key, value in neo4j_id_2_legacy_id_2.items())

        leg_ids_1 = set(neo4j_id_2_legacy_id_1[key] for key in list(neo4j_id_2_matrix_index_1.keys()))
        leg_ids_2 = set(neo4j_id_2_legacy_id_2[key] for key in list(neo4j_id_2_matrix_index_2.keys()))

        first_but_not_second = leg_ids_1 - leg_ids_2
        second_but_not_first = leg_ids_2 - leg_ids_1

        unfold_dict = {}
        # for node_id in chain(first_but_not_second, second_but_not_first):
        #     print node_id
        #     node = DatabaseGraph.find({'legacyID': node_id})[0]
        #     unfold_dict[node_id] = (node_id, node['displayName'],
        #     node.get('forbidden'))
        #
        # first_but_not_second = [unfold_dict[node] for node in first_but_not_second]
        # second_but_not_first = [unfold_dict[node] for node in second_but_not_first]

        log.info('nodes indexed by first laplacian but not second: %s' % first_but_not_second)
        log.info('nodes indexed by second laplacian but not first: %s' % second_but_not_first)

        first_and_second = leg_ids_1.intersection(leg_ids_2)
        input('press enter to continue')

        for legacy_id in first_and_second:
            idx1 = neo4j_id_2_matrix_index_1[legacy_id_2_neo4j_id_1[legacy_id]]
            idx2 = neo4j_id_2_matrix_index_2[legacy_id_2_neo4j_id_2[legacy_id]]

            connections_1 = set(neo4j_id_2_legacy_id_1[matrix_index_2_neo4j_id_1[idx]] for idx in
                             adjacency_matrix_1[idx1, :].nonzero()[1].tolist())
            connections_2 = set(neo4j_id_2_legacy_id_2[matrix_index_2_neo4j_id_2[idx]] for idx in
                             adjacency_matrix_2[idx2, :].nonzero()[1].tolist())

            cons_f_n_s = connections_1 - connections_2
            cons_s_n_f = connections_2 - connections_1

            cons_f_n_s = [leg_id for leg_id in cons_f_n_s if not DatabaseGraph.check_connection_permutation(legacy_id, leg_id)]

            if len(cons_f_n_s):

                if not legacy_id in list(unfold_dict.keys()):
                    node = DatabaseGraph.find({'legacyID': legacy_id})[0]
                    unfold_dict[legacy_id] = (
                        legacy_id, node['displayName'], node.get('forbidden', False))

                for node_id in chain(cons_f_n_s):
                    if not node_id in list(unfold_dict.keys()):
                        node = DatabaseGraph.find({'legacyID': node_id})[0]
                        unfold_dict[node_id] = (
                        node_id, node['displayName'], node.get('forbidden', False))

                cons_f_n_s = [unfold_dict[node] for node in cons_f_n_s if not unfold_dict[node][1]]
                # cons_s_n_f = [unfold_dict[node] for node in cons_s_n_f if not unfold_dict[node][2]]

                # if len(cons_f_n_s) != 0 or len(cons_s_n_f) != 0:
                #     log.info(
                #         'links for %s in first laplacian: %s in second laplacian: %s' % (legacy_id,
                #                                                                                  connections_1,
                #                                                                                  connections_2))

                if len(cons_f_n_s) != 0:
                    log.info(
                        'for %s, connected proteins in first laplacian but not second: %s' % (unfold_dict[legacy_id],
                        cons_f_n_s))

                # if len(cons_s_n_f) != 0:
                #     log.info(
                #         'for %s, connected proteins in second laplacian but not first: %s' % (unfold_dict[legacy_id],
                #         cons_s_n_f))


if __name__ == "__main__":
    interactome_interface_instance = InteractomeInterface()

    interactome_interface_instance.compare_dumps(
        os.path.join(internal_storage, 'mats_compare/old_mat'),
        os.path.join(internal_storage, 'mats_compare/new_mats_4')
    )

    # background_set.full_rebuild()
    # background_set.fast_load()
    # print background_set.pretty_time()
    # print background_set.md5_hash()

    # background_set.set_Uniprot_source(test_set)
    # background_set.compute_current_and_potentials()
    # background_set.export_conduction_system()
    # background_set.randomly_sample([100,250],[5,5], sparse_rounds=10)
