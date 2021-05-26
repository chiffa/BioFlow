"""
Contains all the tools necessary to map GO ontology and Pathway classification from the
database to an Adjacency and Laplacian graph.
"""
import hashlib
import json
import pickle
import random
import string
import math
import datetime
from collections import defaultdict
from copy import copy
from random import shuffle, sample
from csv import reader
from itertools import combinations, chain
from pprint import PrettyPrinter
from random import shuffle
from time import time
import traceback, sys
from typing import Union, Tuple, List

import numpy as np
from scipy.sparse import lil_matrix, triu
from scipy.sparse.csgraph import shortest_path

from bioflow.algorithms_bank import conduction_routines as cr
from bioflow.configs import main_configs as confs
from bioflow.configs.main_configs import Dumps, NewOutputs
from bioflow.sample_storage.mongodb import insert_annotome_rand_samp
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.utils.gdfExportInterface import GdfExportInterface
from bioflow.utils.io_routines import dump_object, undump_object, get_background_bulbs_ids
from bioflow.utils.log_behavior import get_logger
from bioflow.algorithms_bank import sampling_policies
from bioflow.algorithms_bank.flow_calculation_methods import general_flow,\
    reduce_and_deduplicate_sample, evaluate_ops, reduce_ops
from bioflow.algorithms_bank.sampling_policies import characterize_flow_parameters, _is_int


log = get_logger(__name__)


# a pair of debug functions
def _characterise(_object):
    print('Object of size %s and type %s' % (len(_object), type(_object)))


def _characterise_mat(matrix):
    print('Matrix of shape %s, type %s and has %s non-zero terms, min is %s, max is %s' % \
          (matrix.shape, type(matrix), len(matrix.nonzero()[0]), '<>', '<>'))


class GeneOntologyInterface(object):
    """
    Interface between annotome in the knowledge database and the annotome graph laplacian. It is
    heavily skewed towards the Gene Ontology, although can be adapted to be more general than that.

    :param namespace_filter: which namespaces will be used from the annotome (by default the
        "biological process" of the Gene Ontology)
    :param background: (optional) background that will be used for sampling of random
        nodes to build a comparison interface for the
    :param correction_factor:informativity of the node computation correction factors
        (information entropy-wise). (Multiplicative correction factor, additive correction factor)
    :param ultraspec_clean: if the terms considred too specific are excluded
    :param ultraspec_lvl: how many uniprots have to be annotated by a term (directly or
        indirectly) for it not to be considered too specific
    """
    # REFACTOR: [BKI normalization]: move to neo4j parse/insertion types
    _go_up_types = ["is_a_go", "is_part_of_go"]
    _go_reg_types = ["is_Regulant"]

    def __init__(self,
                 namespace_filter=confs.env_bki_filter,
                 background=(),
                 correction_factor=confs.env_bki_correlation_factors,
                 ultraspec_clean=confs.env_bki_ultraspec_clean,
                 ultraspec_lvl=confs.env_bki_ultraspec_lvl):

        self.go_namespace_filter = list(namespace_filter)
        self._background = background
        log.debug('_background set to %d' % len(background))
        self.correction_factor = correction_factor
        self.ultraspec_cleaned = ultraspec_clean
        self.ultraspec_lvl = ultraspec_lvl
        self.init_time = time()
        self.partial_time = time()

        self.entity_2_terms_neo4j_ids = defaultdict(list)
        self.known_up_ids = set()
        self.term_2_entities_neo4j_ids = defaultdict(list)
        self.all_nodes_neo4j_ids = []
        self.node_id_2_mat_idx = {}
        self.mat_idx_2_note_id = {}
        self.total_entropy = None

        self._limiter_reachable_nodes_dict = {}

        self._limiter_up_2_go_reachable_nodes = {}
        self._limiter_go_2_up_reachable_nodes = {}
        self._limiter_up_2_go_step_reachable_nodes = {}
        self._limiter_go_2_up_step_reachable_nodes = {}
        self._limiter_go_2_weighted_ent = {}

        self.neo4j_id_2_display_name = {}
        self.neo4j_id_2_legacy_id = {}
        self.legacy_id_2_neo4j_id = {}

        self.up_neo4j_id_2_leg_id_disp_name = {}

        self.adjacency_matrix = np.zeros((2, 2))
        self.dir_adj_matrix = np.zeros((2, 2))
        self.laplacian_matrix = np.zeros((2, 2))

        self.inflated_laplacian = np.zeros((2, 2))
        self.inflated_idx2lbl = {}
        self.inflated_lbl2idx = {}
        self.binding_intensity = 0

         # REFACTOR [stateless]: this needs to be passed as an argument, not a persistent variable
        self._active_up_sample: List[int] = []
        self._active_weighted_sample: List[Tuple[int, float]] = []
        self._secondary_weighted_sample: Union[None, List[Tuple[int, float]]] = None

        self._flow_calculation_method = general_flow
        self._ops_evaluation_method = evaluate_ops
        self._ops_reduction_method = reduce_ops


        self.UP2UP_voltages = {}
        self.uniprots_2_voltage = {}

        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}

        char_set = string.ascii_uppercase + string.digits
        self.thread_hex = ''.join(random.sample(char_set * 6, 6))

        self.indep_lapl = np.zeros((2, 2))

        log.info('Setting up GO Interface with namespaces %s and %s background UPs',
                 self.go_namespace_filter, len(self._background))

    def pretty_time(self):
        """
        Times the execution

        :return: tuple containing the time since the creation of the Matrix_getter object and
        since the last cal of function formatted as string
        :rtype: str
        """
        it, pt = (round(time() - self.init_time),
                  round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (
            int(it) / 60, it % 60, int(pt) / 60, pt % 60)
        self.partial_time = time()
        return pload

    def _time(self):
        pt = time() - self.partial_time
        return pt

    def reset_thread_hex(self):
        """
        Reset the hex identifier of the object. used when multiprocessing

        :return:
        """
        char_set = string.ascii_uppercase + string.digits
        self.thread_hex = ''.join(sample(char_set * 6, 6))

    def _dump_statics(self):
        dump_object(
            confs.Dumps.GO_builder_stat,
            (self.go_namespace_filter,
             self._background,
             # it does dump the _background from which it will attempt to rebuild itself.
             self.correction_factor,
             self.ultraspec_cleaned,
             self.ultraspec_lvl))

    @staticmethod
    def _undump_statics():
        return undump_object(confs.Dumps.GO_builder_stat)

    def dump_core(self):
        # print(type(self.entity_2_terms_neo4j_ids))
        # print(type(self.term_2_entities_neo4j_ids))
        # print(type(self.deprecated_SeedSet))
        # print(type(self._limiter_reachable_nodes_dict))
        # print(type(self.neo4j_id_2_display_name))
        # print(type(self.neo4j_id_2_legacy_id))
        # print(type(self.legacy_id_2_neo4j_id))
        # print(type(self.all_nodes_neo4j_ids))
        # print(type(self.node_id_2_mat_idx))
        # print(type(self.mat_idx_2_note_id))
        # print(type(self.up_neo4j_id_2_leg_id_disp_name))
        # print(type(self.deprecated_UPs_without_GO))

        dump_object(
            confs.Dumps.GO_dump,
            (self.entity_2_terms_neo4j_ids,
             self.term_2_entities_neo4j_ids,
             self._limiter_reachable_nodes_dict,
             self.neo4j_id_2_display_name,
             self.neo4j_id_2_legacy_id,
             self.legacy_id_2_neo4j_id,
             self.all_nodes_neo4j_ids,
             self.node_id_2_mat_idx,
             self.mat_idx_2_note_id,
             self.up_neo4j_id_2_leg_id_disp_name))

    def _undump_core(self):
        self.entity_2_terms_neo4j_ids, self.term_2_entities_neo4j_ids, self._limiter_reachable_nodes_dict, \
        self.neo4j_id_2_display_name, self.neo4j_id_2_legacy_id, self.legacy_id_2_neo4j_id, self.all_nodes_neo4j_ids, \
        self.node_id_2_mat_idx, self.mat_idx_2_note_id, self.up_neo4j_id_2_leg_id_disp_name =\
            undump_object(confs.Dumps.GO_dump)

        self.known_up_ids = self.entity_2_terms_neo4j_ids.keys()

    def _dump_matrices(self):
        dump_object(
            confs.Dumps.GO_Mats,
            (self.adjacency_matrix,
             self.dir_adj_matrix,
             self.laplacian_matrix))

    def _undump_matrices(self):
        self.adjacency_matrix, self.dir_adj_matrix, self.laplacian_matrix = undump_object(
            confs.Dumps.GO_Mats)

    def _dump_informativities(self):
        dump_object(
            confs.Dumps.GO_Infos,
            (self._limiter_up_2_go_reachable_nodes,
             self._limiter_go_2_up_reachable_nodes,
             self._limiter_up_2_go_step_reachable_nodes,
             self._limiter_go_2_up_step_reachable_nodes,
             self.GO2_Pure_Inf,
             self._limiter_go_2_weighted_ent))

    def _undump_informativities(self):
        self._limiter_up_2_go_reachable_nodes, self._limiter_go_2_up_reachable_nodes, self._limiter_up_2_go_step_reachable_nodes, \
        self._limiter_go_2_up_step_reachable_nodes, self.GO2_Pure_Inf, self._limiter_go_2_weighted_ent = \
            undump_object(confs.Dumps.GO_Infos)

    def _dump_inflated_elements(self):
        dump_object(
            confs.Dumps.GO_Inflated,
            (self.inflated_laplacian,
             self.inflated_idx2lbl,
             self.inflated_lbl2idx,
             self.binding_intensity))

    def _undump_inflated_elements(self):
        self.inflated_laplacian, self.inflated_idx2lbl, \
        self.inflated_lbl2idx, self.binding_intensity = \
            undump_object(confs.Dumps.GO_Inflated)

    def _dump_memoized(self):
        md5 = hashlib.md5(
            json.dumps(sorted(self._active_up_sample), sort_keys=True).encode('utf-8')).hexdigest()

        payload = {
            'UP_hash': md5,
            'sys_hash': self.md5_hash(),
            'UPs': pickle.dumps(self._active_up_sample),
            'currents': pickle.dumps((self.current_accumulator, self.node_current)),
            'voltages': pickle.dumps(self.uniprots_2_voltage)}
        dump_object(confs.Dumps.GO_Analysis_memoized, payload)

    @staticmethod
    def _undump_memoized():
        """
        :return: undumped memoized analysis
        """
        return undump_object(confs.Dumps.GO_Analysis_memoized)

    def _dump_independent_linear_sets(self):
        dump_object(confs.Dumps.GO_Indep_Linset, self.indep_lapl)

    def _undump_independent_linear_sets(self):
        self.indep_lapl = undump_object(confs.Dumps.GO_Indep_Linset)

    def full_rebuild(self):
        """
        Performs a complete rebuild of the InterfaceClass Instance based on parameters provided
        upon construction based on the data in the knowledge database. Upon rebuild saves a copy
        that can be rapidly resurrected with the fast_load() method

        :return: None
        """
        self.annotome_access_and_structure()
        self.get_go_adjacency_and_laplacian()
        self.get_go_reach()
        if self.ultraspec_cleaned:
            self.filter_out_too_specific()
        self.get_laplacians()
        self.inflate_matrix_and_indexes()

        self._dump_statics()
        self.dump_core()
        self._dump_matrices()
        self._dump_informativities()
        self._dump_inflated_elements()

        if self._background:
            if _is_int(self._background[0]):
                self._background = list(set(self.known_up_ids).intersection(set(self._background)))
            else:
                self._background = [(_id, _weight)
                                    for _id, _weight in self._background
                                    if _id in self.known_up_ids]

        else:
            self._background = list(self.known_up_ids)

        log.info('Finished rebuilding the GO Interface object %s', self.pretty_time())

    def fast_load(self):
        """
        Rapidly resurrects the InterfaceClass Instance based on parameters provided
        upon construction. If parameters are mismatched, raises exceptions signalling what
        parameters were mismatched., Trims the background provided upon construction down to what
        actually be sampled (self._background)

        :raise Exception:  wrong filtering namespace parameter
        :raise Exception:  wrong correction factor parameter
        :raise Exception:  wrong ultraspec cleaned parameter
        :raise Exception:  wrong ultraspec level parameter

        """
        namespace_filter, initial_set, correction_factor, ultraspec_cleaned, ultraspec_lvl = \
            self._undump_statics()
        if self.go_namespace_filter != namespace_filter:
            log.critical("Wrong Filtering attempted to be recovered from storage.\n"
                         "\tsaved: %s\n"
                         "\tcurrently active: %s" % (namespace_filter, self.go_namespace_filter))
            raise Exception(
                "Wrong Filtering attempted to be recovered from storage")
        if self.correction_factor != correction_factor:
            log.critical("Wrong correction factor attempted to be recovered from storage")
            raise Exception(
                "Wrong correction factor attempted to be recovered from storage")
        if self.ultraspec_cleaned != ultraspec_cleaned:
            log.critical(
                "Ultraspecific terms leveling state is not the same in the database as requested")
            raise Exception(
                "Ultraspecific terms leveling state is not the same in the database as requested")
        if self.ultraspec_lvl != ultraspec_lvl:
            log.critical(
                "Ultraspecific terms leveling cut-off is not the same in the database as requested")
            raise Exception(
                "Ultraspecific terms leveling cut-off is not the same in the database as requested")

        self._undump_core()

        log.info("_background: %d, entity_2_terms_neo4j_ids %s" % (len(self._background),
                                                                   len(self.known_up_ids)))

        if self._background:
            if _is_int(self._background[0]):
                self._background = list(set(self.known_up_ids).intersection(
                    set(self._background)))
            else:
                self._background = [(_id, _weight)
                                    for _id, _weight in self._background
                                    if _id in self.known_up_ids]

        else:
            self._background = list(self.known_up_ids)

        self._undump_matrices()
        self._undump_informativities()
        self._undump_inflated_elements()


    def annotome_access_and_structure(self, ontology_source=('Gene Ontology')):
        """
        Loads the relationship betweenm the UNIPROTS and annotome as one giant dictionary,
        then between the GO terms themselves

        :param ontology_source:
        :return:
        """
        all_nodes_dict, edges_list = DatabaseGraph.parse_knowledge_entity_net()
        term_counter = 0

        self._limiter_reachable_nodes_dict = defaultdict(lambda: (set(), set(), set(), set()))
        # basically (pure up, out reg, pure down, in_reg)

        for node_id, node_obj in all_nodes_dict.items():
            # uniprot parse
            if list(node_obj.labels)[0] == 'UNIPROT':
                self.up_neo4j_id_2_leg_id_disp_name[node_id] = [node_obj['legacyID'],
                                                                node_obj['displayName']]
            # ontology parse
            else:
                if ontology_source \
                        and node_obj['source'] not in ontology_source:
                    continue
                if self.go_namespace_filter \
                        and node_obj['Namespace'] not in self.go_namespace_filter:
                    continue

                self.neo4j_id_2_display_name[node_id] = node_obj['displayName']
                self.neo4j_id_2_legacy_id[node_id] = node_obj['legacyID']
                self.legacy_id_2_neo4j_id[node_obj['legacyID']] = node_id

                self.all_nodes_neo4j_ids.append(node_id)  # there are also nodes that are annotated by GO
                self.mat_idx_2_note_id[term_counter] = node_id
                self.node_id_2_mat_idx[node_id] = term_counter

                term_counter += 1

        for rel_obj in edges_list:

            start_id = rel_obj.start_node.id
            end_id = rel_obj.end_node.id
            # link uniprots with GO annotations
            # except the linkage does not matter, because it can be linked by other sources
            if ontology_source and all_nodes_dict[end_id]['source'] not in ontology_source:
                continue

            if ontology_source and all_nodes_dict[start_id]['parse_type'] == 'annotation' \
                and all_nodes_dict[start_id]['source'] not in ontology_source:
                continue

            # Uniprot will always be first.
            if self.go_namespace_filter \
                    and all_nodes_dict[end_id]['Namespace'] not in self.go_namespace_filter:
                continue

            # however, both annotations need to be of correct namespace.
            if self.go_namespace_filter \
                    and all_nodes_dict[start_id]['parse_type'] == 'annotation'\
                    and all_nodes_dict[start_id]['Namespace'] not in self.go_namespace_filter:
                continue

            if rel_obj['parse_type'] == 'annotates':
                self.term_2_entities_neo4j_ids[end_id].append(start_id)  # because uniprots are first
                self.entity_2_terms_neo4j_ids[start_id].append(end_id)

            # link annotations between them:
            else:  # el_obj['parse_type'] == 'annotation_relationship':
                # for that we actually will need to build a more complicated
                # OPTIMIZE: to match the previous way it functioned we would have needed to
                #  more up only, not down/sideways. Basically find all the UPs and run the cycle
                #  of rel>GO>rel>GO>rel>GO maps until we are out of Uniprots.
                # OPTIMIZE: that would also allow us to eliminate the overly complex
                #  self.get_go_reach
                #  That would define:
                #   - self._limiter_go_2_up_reachable_nodes
                #   - self._limiter_go_2_up_step_reachable_nodes
                #   - self._limiter_up_2_go_reachable_nodes
                #   - self._limiter_up_2_go_step_reachable_nodes
                #   - GO2_Pure_Inf
                #   - _limiter_go_2_weighted_ent
                # The final decision is that to save the time we will stick with what
                #  there was already before.
                if rel_obj.type in self._go_up_types:
                    self._limiter_reachable_nodes_dict[start_id][0].add(end_id)
                    self._limiter_reachable_nodes_dict[end_id][2].add(start_id)
                elif rel_obj.type in self._go_reg_types:
                    self._limiter_reachable_nodes_dict[start_id][1].add(end_id)
                    self._limiter_reachable_nodes_dict[end_id][3].add(start_id)


        self._limiter_reachable_nodes_dict = dict(self._limiter_reachable_nodes_dict)

        for key in self._limiter_reachable_nodes_dict.keys():
            payload = self._limiter_reachable_nodes_dict[key]
            self._limiter_reachable_nodes_dict[key] = (list(payload[0]),
                                                       list(payload[1]),
                                                       list(payload[2]),
                                                       list(payload[3]))

        self.term_2_entities_neo4j_ids = dict(self.term_2_entities_neo4j_ids)
        self.entity_2_terms_neo4j_ids = dict(self.entity_2_terms_neo4j_ids)

        self.known_up_ids = self.entity_2_terms_neo4j_ids.keys()

    def get_go_adjacency_and_laplacian(self, include_reg=True):
        """
        Builds Undirected and directed adjacency matrices for the GO set and

        :param include_reg: if True, the regulation set is included into the matrix
        :warning: if the parameter above is set to False, get_GO_reach module will be
        unable to function.
        """

        def build_adjacency():
            """
            Builds undirected adjacency matrix for the GO transitions

            """
            base_matrix = lil_matrix((len(self.all_nodes_neo4j_ids), len(self.all_nodes_neo4j_ids)))
            # REFACTOR: [BKI normalization] "package" should be a named tuple
            for node, package in self._limiter_reachable_nodes_dict.items():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.node_id_2_mat_idx[node], self.node_id_2_mat_idx[node2])
                    base_matrix[idx] = 1
                    idx = (idx[1], idx[0])
                    base_matrix[idx] = 1

            self.adjacency_matrix = copy(base_matrix)

        def build_dir_adj():
            """
            Builds directed adjacency matrix for the GO transitions

            """
            base_matrix = lil_matrix((len(self.all_nodes_neo4j_ids), len(self.all_nodes_neo4j_ids)))
            for node, package in self._limiter_reachable_nodes_dict.items():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.node_id_2_mat_idx[node], self.node_id_2_mat_idx[node2])
                    base_matrix[idx] = 1

            self.dir_adj_matrix = copy(base_matrix)

        build_adjacency()
        build_dir_adj()

    def calculate_informativity(self, number, key=None):
        """
        returns an entropy given by a number of equi-probable events, where event is the number.

        :param number:
        """
        # REFACTOR: [Better weights]: rewrite informativity/weight calculation according to a
        #  loadable policy function

        if not self.total_entropy:
            self.total_entropy = - \
                math.log(1. / len(self.known_up_ids), 2)

        if number < 1.0:
            # It actually possible now, the results just won't be used anymore
            log.debug("Term (%s) without reach (%.2f) encountered in informativity calculation "
                      % (key, number))
            return 10 * self.total_entropy
            # raise Exception("Wrong value (%.2f) provided for entropy computation of %s" % (number,
            #                                                                                key)

        if number == 1.0:
            return 2 * self.total_entropy

        return pow(-self.correction_factor[0] * self.total_entropy /
                   math.log(1 / float(number), 2), self.correction_factor[1])


    # REFACTOR: [Maintenability]: method is excessively complex (cyc. complexity ~ 18).
    def get_go_reach(self):
        """
        Recovers by how many different uniprots each GO term is reached, both in
        distance-agnostic and distance-specific terms.

        :raise Exception: if the reaches are not equivalent
        :raise Exception: if null section in a non-nul patch
        """

        def verify_equivalence_of_reaches(step_reach, reach):
            """

            :param step_reach:
            :param reach:
            :raise Exception: if the reaches are not equivalent
            """
            dict_len = {key: [len(val), len(list(step_reach[key].keys()))]
                        for key, val in reach.items()}

            for key, val in dict_len.items():
                if val[1] != val[0]:
                    log.critical(
                        'Reach exploration results not equivalent! Please report the error.')
                    raise Exception(
                        'Reach exploration results not equivalent! Please report the error.')

        def special_sum(_val_dict, filter_function=lambda x: x + 1.0):
            """
            Special sum used for the computation of staged informativity of different terms

            :param _val_dict:
            :param filter_function:
            :raise Exception:
            """
            summer = 0
            for key, val_list in _val_dict.items():
                summer += filter_function(key) * len(val_list)
            return summer

        dir_reg_path = shortest_path(self.dir_adj_matrix, directed=True, method='D')
        dir_reg_path[np.isinf(dir_reg_path)] = 0.0   # potential problem from pycharm
        dir_reg_path = lil_matrix(dir_reg_path)

        # Load all the GOs that can potentially be reached
        self._limiter_go_2_up_reachable_nodes = dict((el, []) for el in list(self._limiter_reachable_nodes_dict.keys()))
        # Load all the UPs that are reached directly from the GO nodes.
        self._limiter_go_2_up_reachable_nodes.update(self.term_2_entities_neo4j_ids)

        pre_go2up_step_reachable_nodes = dict((key, dict((v, 0) for v in val))
                                                for key, val in self._limiter_go_2_up_reachable_nodes.items())
        # when called on possibly un-encoutenred items, anticipate a default
        # value of 100 000

        # Now just scan vertical columns and add UP terms attached
        for idx1, idx2 in zip(list(dir_reg_path.nonzero()[0]),
                              list(dir_reg_path.nonzero()[1])):
            # add UPs annotated by a node to all more general terms.
            self._limiter_go_2_up_reachable_nodes[self.mat_idx_2_note_id[idx2]] +=\
                self._limiter_go_2_up_reachable_nodes[self.mat_idx_2_note_id[idx1]]

            if dir_reg_path[idx1, idx2] < 1.0:
                log.critical("null in non-null patch")
                raise Exception("null in non-null patch")

            step_reach_upgrade = dict(
                (key, val + dir_reg_path[idx1, idx2])
                for key, val in pre_go2up_step_reachable_nodes[self.mat_idx_2_note_id[idx1]].items())

            for k, v in step_reach_upgrade.items():
                pre_go2up_step_reachable_nodes[
                    self.mat_idx_2_note_id[idx2]][k] = min(
                    pre_go2up_step_reachable_nodes[
                        self.mat_idx_2_note_id[idx2]].setdefault(
                        k, 100000), v)

        for key, val in self._limiter_go_2_up_reachable_nodes.items():
            self._limiter_go_2_up_reachable_nodes[key] = list(set(val))

        verify_equivalence_of_reaches(
            pre_go2up_step_reachable_nodes,
            self._limiter_go_2_up_reachable_nodes)

        # Now we need to invert the reach to get the set of all the primary and
        # derived GO terms that describe a UP
        self._limiter_up_2_go_reachable_nodes = dict(
            (key, []) for key in self.known_up_ids)
        self._limiter_up_2_go_step_reachable_nodes = dict(
            (key, defaultdict(list)) for key in self.known_up_ids)
        self._limiter_go_2_up_step_reachable_nodes = dict(
            (key, defaultdict(list)) for key in list(pre_go2up_step_reachable_nodes.keys()))

        for key, val_dict in pre_go2up_step_reachable_nodes.items():
            for k, v in val_dict.items():
                self._limiter_go_2_up_step_reachable_nodes[key][v].append(k)
                self._limiter_up_2_go_step_reachable_nodes[k][v].append(key)
                self._limiter_up_2_go_reachable_nodes[k].append(key)

        # and finally we compute the pure and weighted informativity for each
        # term
        self.GO2_Pure_Inf = dict((key, self.calculate_informativity(len(val), key))
                                 for key, val in self._limiter_go_2_up_reachable_nodes.items())
        self._limiter_go_2_weighted_ent = dict((key, self.calculate_informativity(special_sum(val_dict)))
                                               for key, val_dict in self._limiter_go_2_up_step_reachable_nodes.items())

    def get_laplacians(self):
        """
        Recovers the Laplacian (information conductance) matrixes for the GO annotation terms.
        For weighted laplacian, currently implements a Max-Ent with custom factor as transition
        price.

        :warning: for this method to function, get_GO reach function must be run first.
        :warning: accounting for regulatory relation relation between the GO terms is performed
        if has been done in the adjunction matrix computation

        """
        base_matrix = -copy(self.dir_adj_matrix)
        nz_list = copy(
            list(zip(list(base_matrix.nonzero()[0]), list(base_matrix.nonzero()[1]))))

        # REFACTOR: [Better weights]: change that to a version using a function to calculate the
        #  weights (weigting policy)
        for idx1, idx2 in nz_list:
            min_inf = min(
                self.GO2_Pure_Inf[self.mat_idx_2_note_id[idx1]],
                self.GO2_Pure_Inf[self.mat_idx_2_note_id[idx2]])
            base_matrix[idx1, idx2] = -min_inf
            base_matrix[idx2, idx1] = -min_inf
            base_matrix[idx2, idx2] += min_inf
            base_matrix[idx1, idx1] += min_inf

        self.laplacian_matrix = base_matrix

    def compute_uniprot_dict(self):
        """
        Unused.

        Computes the uniprot method required by a third-party method

        :return:
        """
        uniprot_dict = {}

        for elt in self.up_neo4j_id_2_leg_id_disp_name.keys():
            node = DatabaseGraph.get(elt, 'UNIPROT')
            alt_id = node['legacyID']
            uniprot_dict[alt_id] = (
                elt, self.up_neo4j_id_2_leg_id_disp_name[elt][1])
            uniprot_dict[elt] = alt_id
        pickle.dump(uniprot_dict, open(confs.Dumps.Up_dict_dump, 'wb'))

        return uniprot_dict

    def filter_out_too_specific(self):
        """
        Filters out GO terms that are too specific and builds a directed, undirected adjacency
         maps and laplacian.

        """
        rep_val = self.calculate_informativity(self.ultraspec_lvl)
        self.ultraspec_cleaned = True
        ultraspec_go_terms = list(GO
                                  for GO, reach
                                  in self._limiter_go_2_up_reachable_nodes.items()
                                  if len(reach) < self.ultraspec_lvl)
        for GO in ultraspec_go_terms:
            self.GO2_Pure_Inf[GO] = rep_val

    def md5_hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been
        correctly defined before dump/retrieval
        """
        sorted_initset = sorted(self.node_id_2_mat_idx.keys())

        data = [
            sorted_initset,
            self.go_namespace_filter,
            self.correction_factor,
            self.ultraspec_cleaned,
            self.ultraspec_lvl,
            confs.line_loss,
            confs.use_normalized_laplacian,
            confs.fraction_edges_dropped_in_laplacian]

        md5 = hashlib.md5(json.dumps(data, sort_keys=True).encode('utf-8')).hexdigest()

        log.debug("System md5 hashing done. hash: %s. parameters: \n"
                  "\t sorted init set hash: %s \n"
                  "\t bki environment paramters: %s|%s|%s|%s\n"
                  "\t line loss: %s\n"
                  "\t l_norm: %s\n"
                  "\t edge drop: %s\n"
                  "\t skipping reactome|hint|biogrid: %s\n"
                  % (md5,
                     hashlib.md5(json.dumps(sorted_initset, sort_keys=True).encode('utf-8')).hexdigest(),
                     self.go_namespace_filter,
                     self.correction_factor,
                     self.ultraspec_cleaned,
                     self.ultraspec_lvl,
                     confs.line_loss,
                     confs.use_normalized_laplacian,
                     confs.fraction_edges_dropped_in_laplacian,
                     (confs.env_skip_reactome, confs.env_skip_hint, confs.env_skip_biogrid)))

        return str(md5)

    def active_sample_md5_hash(self, sparse_rounds):
        """
        Performs a hash of characteristics of loaded primary hits list, secondary hits list,
        and background with flow calculation methods. Basically, everything needed to know if a
        random sample is relevant to the currently loaded sample

        :param sparse_rounds: -1 if dense flow calculation, otherwise sparse sampling parameter
        :return:
        """

        sys_hash = self.md5_hash()

        background = []

        if self._background:
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
        """
        Sets the sample to analyze - primary and secondary sources

        :param sample: primary sample being loaded
        :param secondary_sample: secondary sample being loaded (None if none)
        :return:
        """

        def _verify_uniprot_ids(id_weight_vector: List[Tuple[int, float]]):
            uniprots = np.array(id_weight_vector)[:, 0].astype(np.int).tolist()

            if not set(uniprots) <= self.known_up_ids:

                log.warn('Following reached uniprots neo4j_ids were not retrieved upon the '
                         'circulation matrix construction: \n %s',
                         (set(uniprots) - self.known_up_ids))

            _filter = [True
                       if uniprot in self.known_up_ids
                       else False
                       for uniprot in uniprots]

            return np.array(id_weight_vector)[_filter, :].tolist()

        self._active_weighted_sample = _verify_uniprot_ids(reduce_and_deduplicate_sample(sample))

        self._active_up_sample = np.array(self._active_weighted_sample)[:, 0].tolist()

        if secondary_sample is not None:

            log.debug('debug: secondary_weight sample %s' % secondary_sample)

            self._secondary_weighted_sample = \
                _verify_uniprot_ids(reduce_and_deduplicate_sample(secondary_sample))

            # KNOWNBUG: logic bug occurs if the ids supplied have no annotations in the db

            log.debug('debug: secondary_weight sample %s' % np.array(
                self._secondary_weighted_sample))

            self._active_up_sample = list(set(self._active_up_sample
                                              + np.array(self._secondary_weighted_sample)[:, 0].tolist()))

    def evaluate_ops(self, sparse_rounds=-1):
        """
        Evaluate the number of pairs of nodes that wlll be used for flow calculation

        :param sparse_rounds: sparse rounds parameter, if -1 will be considered dense
        :return:
        """
        log.debug('evaluate_ops call')
        ro = sampling_policies.characterize_flow_parameters(self._active_weighted_sample,
                                                            self._secondary_weighted_sample,
                                                            -2)
        return self._ops_evaluation_method(ro[1], ro[3], sparse_rounds)

    def reduce_ops(self, ops_limit):
        """
        Evaluates the value of the sparse_round parameter need to keep the number of pairs of
        nodes used for flow calculation under a given budget

        :param ops_limit: node pair budget
        :return:
        """
        log.debug('reduce_ops call')
        ro = sampling_policies.characterize_flow_parameters(self._active_weighted_sample,
                                                            self._secondary_weighted_sample,
                                                            -2)

        return self._ops_reduction_method(ro[1], ro[3], ops_limit)

    def inflate_matrix_and_indexes(self):
        """
        Performs the laplacian matrix inflation to incorporate the uniprots on which we
         will be running the analysis
        """
        # branching distribution: at least 10x the biggest conductivity of the
        # system, unless too specific, in which case ~ specific level
        self.binding_intensity = 10 * self.calculate_informativity(self.ultraspec_lvl)
        fixed_index = self.laplacian_matrix.shape[0]

        up2idxs = dict((UP, fixed_index + Idx)
                       for Idx, UP in enumerate(self.known_up_ids))
        idx2ups = dict((Idx, UP) for UP, Idx in up2idxs.items())

        self.inflated_laplacian = lil_matrix(
            (self.laplacian_matrix.shape[0] + len(self.known_up_ids),
             self.laplacian_matrix.shape[1] + len(self.known_up_ids)))

        self.inflated_laplacian[:self.laplacian_matrix.shape[0], :self.laplacian_matrix.shape[1]] =\
            self.laplacian_matrix

        for uniprot in self.known_up_ids:
            for go_term in self.entity_2_terms_neo4j_ids.get(uniprot, []):  # should never hit the [] though
                self.inflated_laplacian[
                    up2idxs[uniprot], up2idxs[uniprot]] += self.binding_intensity
                self.inflated_laplacian[
                    self.node_id_2_mat_idx[go_term], self.node_id_2_mat_idx[go_term]] += self.binding_intensity
                self.inflated_laplacian[
                    self.node_id_2_mat_idx[go_term], up2idxs[uniprot]] -= self.binding_intensity
                self.inflated_laplacian[
                    up2idxs[uniprot], self.node_id_2_mat_idx[go_term]] -= self.binding_intensity

        self.inflated_lbl2idx = copy(self.node_id_2_mat_idx)
        self.inflated_lbl2idx.update(up2idxs)
        self.inflated_idx2lbl = copy(self.mat_idx_2_note_id)
        self.inflated_idx2lbl.update(idx2ups)

    def set_uniprot_source(self, uniprots):
        """
        Sets the deprecated_reached_uniprots_neo4j_id_list on which the circulation computation routines
        will be performed by the otehr methods.Avoids passing as argument large lists of parameters.

        :param uniprots: List of node IDs of the uniprots on which we would like to perform
        current computations
        :raise Warning: if the uniprots were not present in the set of GOs for which we
        built the system or had no GO attached to them
        """
        if not set(uniprots) <= self.known_up_ids :
            na_set = set(uniprots) - self.known_up_ids
            log.warning('%s uniprots out of %s either were not present in the constructions set '
                        'or have no GO terms attached to them.', len(na_set), len(set(uniprots)))
            log.debug('full list of uniprots that cannot be analyzed: \n%s', na_set)

        self._active_up_sample = [uniprot for uniprot in uniprots if uniprot in self.known_up_ids]


    def compute_current_and_potentials(
            self,
            memoized: bool = True,
            incremental: bool = False,  # This should always be false and was used in order to
            # resume the sampling
            cancellation: bool = False,
            sparse_rounds: int = -1,
            potential_dominated: bool = True):
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier
        knowledge flow analysis

        :param memoized: if the tensions and individual relation matrices should be stored in
         the matrix and dumped at the end computation (required for submatrix re-computation)
        :param incremental: if True, all the circulation computation will be added to the
        existing ones. Useful for the computation of particularly big systems with
        intermediate dumps
        :param cancellation: divides the final current by #Nodes**2/2, i.e. makes the currents
        comparable between circulation systems of different  sizes.
        :param sparse_rounds: if set to a positive integer the sampling will be sparse and
        not dense, i.e. instead of computation for each node pair, only an estimation will be
        made, equal to computing sparse sampling association with other randomly chosen nodes
        :param potential_dominated: if the total current is normalized to potential
        :return: adjusted conduction system
        """

        if not incremental or self.current_accumulator == np.zeros((2, 2)):
            self.current_accumulator = lil_matrix(self.inflated_laplacian.shape)
            self.UP2UP_voltages = {}
            self.uniprots_2_voltage = {}  # REFACTOR [maintenance]: remove


        weighted_up_pairs = self._flow_calculation_method(self._active_weighted_sample,
                                                          self._secondary_weighted_sample,
                                                          sparse_rounds)
        # pairs in the list of pairs are now guaranteed to be weighted

        total_pairs = len(weighted_up_pairs)
        breakpoints = 300
        previous_time = time()

        for counter, (up_w_id_1, up_w_id_2) in enumerate(weighted_up_pairs):

            up_id_1, w_1 = up_w_id_1
            up_id_2, w_2 = up_w_id_2
            mean_weight = (w_1 + w_2) / 2.

            idx1, idx2 = (self.inflated_lbl2idx[up_id_1], self.inflated_lbl2idx[up_id_2])

            pre_reach = self._limiter_up_2_go_reachable_nodes[up_id_1] + \
                        self._limiter_up_2_go_reachable_nodes[up_id_2] + \
                        [up_id_1] + [up_id_2]

            reach = [self.inflated_lbl2idx[label] for label in pre_reach]

            current_upper, potential_diff = cr.group_edge_current_with_limitations(
                inflated_laplacian=self.inflated_laplacian,
                idx_pair=(idx1, idx2),
                reach_limiter=reach)

            self.UP2UP_voltages[tuple(sorted((up_id_1, up_id_2)))] = potential_diff

            if potential_dominated:
                if potential_diff != 0:
                    current_upper = current_upper / potential_diff

                else:
                    log.warning('pairwise flow. On indexes %s %s potential difference is null. %s',
                                up_id_1, up_id_2, 'Tension-normalization was aborted')

            self.current_accumulator = self.current_accumulator + \
                                       cr.sparse_abs(current_upper) * mean_weight

            if counter % breakpoints == 0 and counter > 1:
                # TODO: [load bar] the internal loop load bar goes here
                compops = float(breakpoints) / (time() - previous_time)
                mins_before_termination = (total_pairs - counter) / compops // 60
                finish_time = datetime.datetime.now() + datetime.timedelta(minutes=mins_before_termination)
                log.info("thread hex: %s; progress: %s/%s, current speed: %.2f compop/s, "
                         "time remaining: "
                         "%.0f "
                         "min, finishing: %s "
                         % (self.thread_hex, counter, total_pairs, compops, mins_before_termination,
                            finish_time.strftime("%m/%d/%Y, %H:%M:%S")))
                previous_time = time()

        self.current_accumulator = triu(self.current_accumulator)

        if cancellation:
            self.current_accumulator /= float(total_pairs)

        index_current = cr.get_current_through_nodes(self.current_accumulator)

        log.info('current accumulator shape %s, sum %s',
                 self.current_accumulator.shape, np.sum(self.current_accumulator))

        self.node_current = dict((self.inflated_idx2lbl[idx], val)
                                 for idx, val in enumerate(index_current))

        if memoized:
            self._dump_memoized()

    def format_node_props(self, node_current, limit=0.01):
        """
        Formats the nodes for the analysis by in the knowledge_access_analysis module

        :param node_current: Current through the GO nodes
        :param limit: hard limit to go_namespace_filter out the GO terms with too little current
         (compensates the minor currents in the gird)s
        :return: {GO:[node current, pure GO informativity, Number of reachable nodes]}
        """
        characterization_dict = {}
        limiting_current = max(node_current.values()) * limit
        log.debug('formatting node props with %.2f limiting current' % limiting_current)

        for go_term in self.node_id_2_mat_idx.keys():

            if node_current[go_term] > limiting_current:
                characterization_dict[go_term] = [
                    node_current[go_term],
                    self.GO2_Pure_Inf[go_term],
                    len(self._limiter_go_2_up_reachable_nodes[go_term])]

        # should never occur, unless single node is a massive outlier
        if len(characterization_dict) < 2:
            characterization_dict = self.format_node_props(node_current, limit/100)

        return characterization_dict

    def export_conduction_system(self,
                                 p_value_dict: dict = None,
                                 output_location: str = ''):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format
         and flushes it into a file that can be viewed with Gephi

        :param p_value_dict:
        :raise Warning:
        """
        node_char_names = [
            'Current',
            'Type',
            'Legacy_ID',
            'Names',
            'Pure_informativity',
            'Confusion_potential',
            'p-value',
            'p_p-value',
            'Source_W']

        node_char_types = [
            'DOUBLE',
            'VARCHAR',
            'VARCHAR',
            'VARCHAR',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE']

        if p_value_dict is None:
            p_value_dict = defaultdict(lambda: (np.nan, np.nan, np.nan))

        nan_neg_log10 = lambda x: x if type(x) is str else -np.log10(x)

        char_dict = {}

        for GO in self.node_id_2_mat_idx.keys():
            char_dict[GO] = [str(self.node_current[GO]),
                             'GO', self.neo4j_id_2_legacy_id[GO],
                             self.neo4j_id_2_display_name[GO].replace(',', '-'),
                             str(self.GO2_Pure_Inf[GO]),
                             str(len(self._limiter_go_2_up_reachable_nodes[GO])),
                             str(p_value_dict[int(GO)][0]),
                             str(nan_neg_log10(p_value_dict[int(GO)][0])),
                             '0']

        for UP in self._active_up_sample:

            in_sample_weight = 1
            for node, weight in self._active_weighted_sample:
                if UP == node:
                    in_sample_weight = weight

            if self._secondary_weighted_sample is not None:
                for node, weight in self._secondary_weighted_sample:
                    if UP == node:
                        in_sample_weight = -weight

            char_dict[UP] = [str(self.node_current[UP]),
                             'UP', self.up_neo4j_id_2_leg_id_disp_name[UP][0],
                             str(self.up_neo4j_id_2_leg_id_disp_name[UP][1]).replace(',', '-'),
                             str(self.binding_intensity),
                             '1',
                             '0.1',
                             '1',
                             str(in_sample_weight)]  # TODOC: document the defaults

        if output_location == '':
            output_location = NewOutputs().GO_GDF_output

        gdf_exporter = GdfExportInterface(
            target_fname=output_location,
            field_names=node_char_names,
            field_types=node_char_types,
            node_properties_dict=char_dict,
            min_current=0.01,
            index_2_label=self.inflated_idx2lbl,
            label_2_index=self.inflated_lbl2idx,
            current_matrix=self.current_accumulator)
        # TODO: [Better stats]: twister compared to random sample?
        gdf_exporter.write()

    def randomly_sample(
            self,
            random_samples,
            sparse_rounds: int = -1,
            no_add=False,
            pool_no=None,
            sampling_policy=sampling_policies.matched_sampling,
            optional_sampling_param = 'exact'):
        """
        Randomly samples the set of deprecated_reached_uniprots_neo4j_id_list used to create the model.
        This is the null model creation routine

        :param random_samples: how many times we would like to sample each unirot number
        :param sparse_rounds:  if we want to use sparse sampling
            (useful in case of large uniprot sets),
        :param no_add: if set to True, the result of sampling will not be added to the database
            of samples. Useful if re-running tests with similar parameters several times.
        :param pool_no: explicit sampling pool number (used for reporting/debugging)
        :param sampling_policy: sampling policy used
        :param optional_sampling_param: sampling policy optional argument
        """

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
                                                     self._background,
                                                     random_samples,
                                                     optional_sampling_param):

            self.set_flow_sources(sample, sec_sample)

            sample_chars = characterize_flow_parameters(sample, sec_sample, sparse_rounds)
            sample_hash = sample_chars[-1]

            log.info('Sampling thread: %s, Thread hex: %s; Random sample %d/%d \n'
                     'sampling characteristics: sys_hash: %s, sample_hash: %s, '
                     'target_hash: %s' %
                     (pool_no, self.thread_hex, i, random_samples,
                      self.md5_hash(), sample_hash, super_hash))

            # TODO: [load bar]: the external loop progress bar goes here

            # TODO: [fast resurrection] fast resurrection is impossible (memoized is false,
            #  but pipeline is broken)
            self.compute_current_and_potentials(memoized=False, sparse_rounds=sparse_rounds)

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

            insert_annotome_rand_samp(
                        {
                            'UP_hash': sample_ids_md5,
                            'sys_hash': self.md5_hash(),
                            'active_sample_hash': self.active_sample_md5_hash(sparse_rounds),
                            'target_sample_hash': super_hash,
                            'sampling_policy': sampling_policy.__name__,
                            'sampling_policy_options': optional_sampling_param,
                            'sparse_rounds': sparse_rounds,
                            'UPs': pickle.dumps(self._active_up_sample),
                            'sample': pickle.dumps(self._active_weighted_sample),
                            'sec_sample': pickle.dumps(self._secondary_weighted_sample),
                            'currents': pickle.dumps(
                                (self.current_accumulator,
                                 self.node_current)),
                            'voltages': pickle.dumps(
                                self.UP2UP_voltages)})

        self._active_weighted_sample = preserved_sample
        self._secondary_weighted_sample = preserved_sec_sample

    def get_independent_linear_groups(self):
        """
        Recovers independent linear groups of the GO terms. Independent linear groups are
        those that share a significant amount of reached_uniprots_neo4j_id_list in common
        """
        self.indep_lapl = lil_matrix((len(self.all_nodes_neo4j_ids), len(self.all_nodes_neo4j_ids)))
        for GO_list in self._limiter_up_2_go_reachable_nodes.values():
            for GO1, GO2 in combinations(GO_list, 2):
                idx1, idx2 = (self.node_id_2_mat_idx[GO1], self.node_id_2_mat_idx[GO2])
                self.indep_lapl[idx1, idx2] += -1
                self.indep_lapl[idx2, idx1] += -1
                self.indep_lapl[idx2, idx2] += 1
                self.indep_lapl[idx1, idx1] += 1


if __name__ == '__main__':
    # Creates an instance of MatrixGetter and loads pre-computed values

    go_interface_instance = GeneOntologyInterface(
        background=get_background_bulbs_ids())
    go_interface_instance.full_rebuild()

    # loading takes 1-6 seconds.
    # fill for reach only is done in 2 seconds,
    # tepping takes another 15,
    # inverting + info computation - 1 more second
    # Laplacian building =>
    ##
    # full computation - 3 minutes 18 seconds; save 7 seconds, retrieval - 3
    # seconds

    # go_interface_instance.fast_load()
    # print go_interface_instance.pretty_time()

    # go_interface_instance.get_indep_linear_groups()
    # go_interface_instance.dump_Indep_Linset()

    # go_interface_instance.randomly_sample([10, 25], [5]*2, chromosome_specific=15)

    # go_interface_instance.set_Uniprot_source(experimental)
    # go_interface_instance.compute_current_and_potentials(sparse_rounds=10)
    # go_interface_instance.export_conduction_system()

    # go_interface_instance.deprecated_export_subsystem(experimental, ['186958', '142401', '147798', '164077'])

    # data_array = np.array([log(val) for val in go