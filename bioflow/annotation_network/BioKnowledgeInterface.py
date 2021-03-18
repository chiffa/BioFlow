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
from csv import reader
from itertools import combinations, chain
from pprint import PrettyPrinter
from random import shuffle
from time import time
import traceback, sys

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

log = get_logger(__name__)


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
    :param background_up_ids: (optional) background that will be used for sampling of random
        nodes to build a comparison interface for the
    :param correction_factor:informativity of the node computation correction factors
        (information entropy-wise). (Multiplicative correction factor, additive correction factor)
    :param ultraspec_clean: if the terms considred too specific are excluded
    :param ultraspec_lvl: how many uniprots have to be annotated by a term (directly or
        indirectly) for it not to be considered too specific
    """
    # REFACTOR: [BKI normalization]: move to configs
    _GOUpTypes = ["is_a_go", "is_part_of_go"]
    _GORegTypes = ["is_Regulant"]

    def __init__(self,
                 namespace_filter=('biological_process',),
                 background_up_ids=(),
                 correction_factor=(1, 1),
                 ultraspec_clean=True,
                 ultraspec_lvl=3):

        self.go_namespace_filter = list(namespace_filter)
        self._background = background_up_ids
        log.debug('_background set to %d' % len(background_up_ids))
        self.correction_factor = correction_factor
        self.ultraspec_cleaned = ultraspec_clean
        self.ultraspec_lvl = ultraspec_lvl
        self.init_time = time()
        self.partial_time = time()

        self.UP2GO_Dict = defaultdict(list)
        self.GO2UP = defaultdict(list)
        self.All_GOs = []
        self.GO2Num = {}
        self.Num2GO = {}
        self.total_entropy = None

        self.reachable_nodes_dict = {}

        self.UP2GO_Reachable_nodes = {}
        self.GO2UP_Reachable_nodes = {}
        self.UP2GO_step_Reachable_nodes = {}
        self.GO2UP_step_Reachable_nodes = {}
        self.GO2_Weighted_Ent = {}

        self.GO_names = {}
        self.GO_legacy_ids = {}
        self.rev_GO_ids = {}

        self.UP_names = {}

        self.adjacency_matrix = np.zeros((2, 2))
        self.dir_adj_matrix = np.zeros((2, 2))
        self.laplacian_matrix = np.zeros((2, 2))

        self.inflated_Laplacian = np.zeros((2, 2))
        self.inflated_idx2lbl = {}
        self.inflated_lbl2idx = {}
        self.binding_intensity = 0

         # REFACTOR [stateless]: this needs to be passed as an argument, not a persistent variable
        self.active_up_sample = []
        self.UP2UP_voltages = {}
        self.uniprots_2_voltage = {}

        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}

        char_set = string.ascii_uppercase + string.digits
        self.thread_hex = ''.join(random.sample(char_set * 6, 6))

        self.indep_lapl = np.zeros((2, 2))
        self.sparsely_sampled = False

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
        # print(type(self.UP2GO_Dict))
        # print(type(self.GO2UP))
        # print(type(self.deprecated_SeedSet))
        # print(type(self.reachable_nodes_dict))
        # print(type(self.GO_names))
        # print(type(self.GO_legacy_ids))
        # print(type(self.rev_GO_ids))
        # print(type(self.All_GOs))
        # print(type(self.GO2Num))
        # print(type(self.Num2GO))
        # print(type(self.UP_names))
        # print(type(self.deprecated_UPs_without_GO))

        dump_object(
            confs.Dumps.GO_dump,
            (self.UP2GO_Dict,
             self.GO2UP,
             self.reachable_nodes_dict,
             self.GO_names,
             self.GO_legacy_ids,
             self.rev_GO_ids,
             self.All_GOs,
             self.GO2Num,
             self.Num2GO,
             self.UP_names))

    def _undump_core(self):
        self.UP2GO_Dict, self.GO2UP, self.reachable_nodes_dict, \
        self.GO_names, self.GO_legacy_ids, self.rev_GO_ids, self.All_GOs, \
        self.GO2Num, self.Num2GO, self.UP_names =\
            undump_object(confs.Dumps.GO_dump)

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
            (self.UP2GO_Reachable_nodes,
             self.GO2UP_Reachable_nodes,
             self.UP2GO_step_Reachable_nodes,
             self.GO2UP_step_Reachable_nodes,
             self.GO2_Pure_Inf,
             self.GO2_Weighted_Ent))

    def _undump_informativities(self):
        self.UP2GO_Reachable_nodes, self.GO2UP_Reachable_nodes, self.UP2GO_step_Reachable_nodes, \
            self.GO2UP_step_Reachable_nodes, self.GO2_Pure_Inf, self.GO2_Weighted_Ent = \
            undump_object(confs.Dumps.GO_Infos)

    def _dump_inflated_elements(self):
        dump_object(
            confs.Dumps.GO_Inflated,
            (self.inflated_Laplacian,
             self.inflated_idx2lbl,
             self.inflated_lbl2idx,
             self.binding_intensity))

    def _undump_inflated_elements(self):
        self.inflated_Laplacian, self.inflated_idx2lbl, \
            self.inflated_lbl2idx, self.binding_intensity = \
            undump_object(confs.Dumps.GO_Inflated)

    def _dump_memoized(self):
        md5 = hashlib.md5(
            json.dumps(sorted(self.active_up_sample), sort_keys=True).encode('utf-8')).hexdigest()

        payload = {
            'UP_hash': md5,
            'sys_hash': self.md5_hash(),
            'size': len(self.active_up_sample),
            'UPs': pickle.dumps(self.active_up_sample),
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

        log.info('Finished rebuilding the GO Interface object %s', self.pretty_time())

    def fast_load(self):
        """
        Rapidly resurrects the InterfaceClass Instance based on parameters provided
        upon construction. If parameters are mismatched, raises exceptions signalling what
        parameters were mismatched., Trims the background provided upon construction down to what
        actually be sampled (self._background)

        """
        namespace_filter, initial_set, correction_factor, ultraspec_cleaned, ultraspec_lvl = \
            self._undump_statics()
        if self.go_namespace_filter != namespace_filter:
            log.critical("Wrong Filtering attempted to be recovered from storage")
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

        log.info("_background: %d, UP2GO_Dict %s" % (len(self._background), len(self.UP2GO_Dict)))

        if not self._background:
            self._background = list(self.UP2GO_Dict.keys())
        else:
            self._background = list(set(self._background).intersection(set(self.UP2GO_Dict.keys())))

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
        # self._background here is irrelevant
        # self.deprecated_SeedSet becomes deprecated as well
        # self._GOUpTypes dissapear
        # self._GO_Reg_Types disappear as we;;

        all_nodes_dict, edges_list = DatabaseGraph.parse_knowledge_entity_net()
        term_counter = 0

        self.reachable_nodes_dict = defaultdict(lambda: (set(), set(), set(), set()))
        # basically (pure up, out reg, pure down, in_reg)

        for node_id, node_obj in all_nodes_dict.items():
            # uniprot parse
            if list(node_obj.labels)[0] == 'UNIPROT':
                self.UP_names[node_id] = [node_obj['legacyID'],
                                          node_obj['displayName']]
            # ontology parse
            else:
                if ontology_source \
                        and node_obj['source'] not in ontology_source:
                    continue
                if self.go_namespace_filter \
                        and node_obj['Namespace'] not in self.go_namespace_filter:
                    continue

                self.GO_names[node_id] = node_obj['displayName']
                self.GO_legacy_ids[node_id] = node_obj['legacyID']
                self.rev_GO_ids[node_obj['legacyID']] = node_id

                self.All_GOs.append(node_id)
                self.Num2GO[term_counter] = node_id
                self.GO2Num[node_id] = term_counter

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
                self.GO2UP[end_id].append(start_id)  # because uniprots are first
                self.UP2GO_Dict[start_id].append(end_id)

            # link annotations between them:
            else:  # el_obj['parse_type'] == 'annotation_relationship':
                # for that we actually will need to build a more complicated
                # OPTIMIZE: to match the previous way it functioned we would have needed to
                #  more up only, not down/sideways. Basically find all the UPs and run the cycle
                #  of rel>GO>rel>GO>rel>GO maps until we are out of Uniprots.
                # OPTIMIZE: that would also allow us to eliminate the overly complex
                #  self.get_go_reach
                #  That would define:
                #   - self.GO2UP_Reachable_nodes
                #   - self.GO2UP_step_Reachable_nodes
                #   - self.UP2GO_Reachable_nodes
                #   - self.UP2GO_step_Reachable_nodes
                #   - GO2_Pure_Inf
                #   - GO2_Weighted_Ent
                # The final decision is that to save the time we will stick with what
                #  there was already before.
                if rel_obj.type in self._GOUpTypes:
                    self.reachable_nodes_dict[start_id][0].add(end_id)
                    self.reachable_nodes_dict[end_id][2].add(start_id)
                elif rel_obj.type in self._GORegTypes:
                    self.reachable_nodes_dict[start_id][1].add(end_id)
                    self.reachable_nodes_dict[end_id][3].add(start_id)


        self.reachable_nodes_dict = dict(self.reachable_nodes_dict)

        for key in self.reachable_nodes_dict.keys():
            payload = self.reachable_nodes_dict[key]
            self.reachable_nodes_dict[key] = (list(payload[0]),
                                              list(payload[1]),
                                              list(payload[2]),
                                              list(payload[3]))

        self.GO2UP = dict(self.GO2UP)
        self.UP2GO_Dict = dict(self.UP2GO_Dict)

        self._background = list(self.UP2GO_Dict.keys())

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
            base_matrix = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
            # REFACTOR: [BKI normalization] "package" should be a named tuple
            for node, package in self.reachable_nodes_dict.items():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.GO2Num[node], self.GO2Num[node2])
                    base_matrix[idx] = 1
                    idx = (idx[1], idx[0])
                    base_matrix[idx] = 1

            self.adjacency_matrix = copy(base_matrix)

        def build_dir_adj():
            """
            Builds directed adjacency matrix for the GO transitions

            """
            base_matrix = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
            for node, package in self.reachable_nodes_dict.items():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.GO2Num[node], self.GO2Num[node2])
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
                math.log(1. / len(self.UP2GO_Dict.keys()), 2)

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


    # REFACTOR: [Sanity]: ethod is excessively complex (cyc. complexity ~ 18).
    def get_go_reach(self):
        """
        Recovers by how many different uniprots each GO term is reached, both in
        distance-agnostic and distance-specific terms.
        """

        def verify_equivalence_of_reaches(step_reach, reach):
            """

            :param step_reach:
            :param reach:
            :raise Exception:
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
        self.GO2UP_Reachable_nodes = dict((el, []) for el in list(self.reachable_nodes_dict.keys()))
        # Load all the UPs that are reached directly from the GO nodes.
        self.GO2UP_Reachable_nodes.update(self.GO2UP)

        pre_go2up_step_reachable_nodes = dict((key, dict((v, 0) for v in val))
                                                for key, val in self.GO2UP_Reachable_nodes.items())
        # when called on possibly un-encoutenred items, anticipate a default
        # value of 100 000

        # Now just scan vertical columns and add UP terms attached
        for idx1, idx2 in zip(list(dir_reg_path.nonzero()[0]),
                              list(dir_reg_path.nonzero()[1])):
            # add UPs annotated by a node to all more general terms.
            self.GO2UP_Reachable_nodes[self.Num2GO[idx2]] +=\
                self.GO2UP_Reachable_nodes[self.Num2GO[idx1]]

            if dir_reg_path[idx1, idx2] < 1.0:
                log.critical("null in non-null patch")
                raise Exception("null in non-null patch")

            step_reach_upgrade = dict(
                (key, val + dir_reg_path[idx1, idx2])
                for key, val in pre_go2up_step_reachable_nodes[self.Num2GO[idx1]].items())

            for k, v in step_reach_upgrade.items():
                pre_go2up_step_reachable_nodes[
                    self.Num2GO[idx2]][k] = min(
                    pre_go2up_step_reachable_nodes[
                        self.Num2GO[idx2]].setdefault(
                        k, 100000), v)

        for key, val in self.GO2UP_Reachable_nodes.items():
            self.GO2UP_Reachable_nodes[key] = list(set(val))

        verify_equivalence_of_reaches(
            pre_go2up_step_reachable_nodes,
            self.GO2UP_Reachable_nodes)

        # Now we need to invert the reach to get the set of all the primary and
        # derived GO terms that describe a UP
        self.UP2GO_Reachable_nodes = dict(
            (key, []) for key in list(self.UP2GO_Dict.keys()))
        self.UP2GO_step_Reachable_nodes = dict(
            (key, defaultdict(list)) for key in list(self.UP2GO_Dict.keys()))
        self.GO2UP_step_Reachable_nodes = dict(
            (key, defaultdict(list)) for key in list(pre_go2up_step_reachable_nodes.keys()))

        for key, val_dict in pre_go2up_step_reachable_nodes.items():
            for k, v in val_dict.items():
                self.GO2UP_step_Reachable_nodes[key][v].append(k)
                self.UP2GO_step_Reachable_nodes[k][v].append(key)
                self.UP2GO_Reachable_nodes[k].append(key)

        # and finally we compute the pure and weighted informativity for each
        # term
        self.GO2_Pure_Inf = dict((key, self.calculate_informativity(len(val), key))
                                 for key, val in self.GO2UP_Reachable_nodes.items())
        self.GO2_Weighted_Ent = dict((key, self.calculate_informativity(special_sum(val_dict)))
                                    for key, val_dict in self.GO2UP_step_Reachable_nodes.items())

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
                self.GO2_Pure_Inf[self.Num2GO[idx1]],
                self.GO2_Pure_Inf[self.Num2GO[idx2]])
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

        for elt in self.UP_names.keys():
            node = DatabaseGraph.get(elt, 'UNIPROT')
            alt_id = node['legacyID']
            uniprot_dict[alt_id] = (
                elt, self.UP_names[elt][1])
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
                                  in self.GO2UP_Reachable_nodes.items()
                                  if len(reach) < self.ultraspec_lvl)
        for GO in ultraspec_go_terms:
            self.GO2_Pure_Inf[GO] = rep_val

    def md5_hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been
        correctly defined before dump/retrieval
        """
        sorted_initset = sorted(self._background)

        data = [
            self.go_namespace_filter,
            sorted_initset,
            self.correction_factor,
            self.ultraspec_cleaned,
            self.ultraspec_lvl]

        md5 = hashlib.md5(json.dumps(data, sort_keys=True).encode('utf-8')).hexdigest()

        return str(md5)

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
                       for Idx, UP in enumerate(self.UP2GO_Dict.keys()))
        idx2ups = dict((Idx, UP) for UP, Idx in up2idxs.items())

        self.inflated_Laplacian = lil_matrix(
            (self.laplacian_matrix.shape[0] + len(set(self.UP2GO_Dict.keys())),
             self.laplacian_matrix.shape[1] + len(set(self.UP2GO_Dict.keys()))))

        self.inflated_Laplacian[:self.laplacian_matrix.shape[0], :self.laplacian_matrix.shape[1]] =\
            self.laplacian_matrix

        for uniprot in self.UP2GO_Dict.keys():
            for go_term in self.UP2GO_Dict.get(uniprot, []):  # should never hit the [] though
                self.inflated_Laplacian[
                    up2idxs[uniprot], up2idxs[uniprot]] += self.binding_intensity
                self.inflated_Laplacian[
                    self.GO2Num[go_term], self.GO2Num[go_term]] += self.binding_intensity
                self.inflated_Laplacian[
                    self.GO2Num[go_term], up2idxs[uniprot]] -= self.binding_intensity
                self.inflated_Laplacian[
                    up2idxs[uniprot], self.GO2Num[go_term]] -= self.binding_intensity

        self.inflated_lbl2idx = copy(self.GO2Num)
        self.inflated_lbl2idx.update(up2idxs)
        self.inflated_idx2lbl = copy(self.Num2GO)
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
        if not set(uniprots) <= set(self.UP2GO_Dict.keys()):
            na_set = set(uniprots) - set(self.UP2GO_Dict.keys())
            log.warning('%s uniprots out of %s either were not present in the constructions set '
                        'or have no GO terms attached to them.', len(na_set), len(set(uniprots)))
            log.debug('full list of uniprots that cannot be analyzed: \n%s', na_set)

        self.active_up_sample = [
            uniprot for uniprot in uniprots if uniprot in list(self.UP2GO_Dict.keys())]

        log.info("tried to set up list %d, intersected with %d UP2GO_dict.keys(), ended up with %d"
                 % (len(uniprots), len(self.UP2GO_Dict.keys()), len(self.active_up_sample)))

    def compute_current_and_potentials(
            self,
            memoized: bool = True,
            incremental: bool = False,  # This should always be false and was used in order to
            # resume the sampling
            cancellation: bool = True,
            sparse_samples=False):
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
        :param sparse_samples: if set to an integer the sampling will be sparse and not dense,
         i.e. instead of computation for each node pair, only an estimation will be made, equal to
         computing sparse_samples association with other randomly chosen nodes
        :return: adjusted conduction system
        """
        if not incremental or self.current_accumulator == np.zeros((2, 2)):
            self.current_accumulator = lil_matrix(self.inflated_Laplacian.shape)
            self.UP2UP_voltages = {}
            self.uniprots_2_voltage = {}

        iterator = []
        if sparse_samples:
            for _ in range(0, sparse_samples):
                _sample = copy(self.active_up_sample)
                random.shuffle(_sample)
                iterator += list(zip(_sample[:len(_sample) // 2], _sample[len(_sample) // 2:]))
                self.sparsely_sampled = True
        else:
            iterator = combinations(self.active_up_sample, 2)

        iterator = [item for item in iterator]

        total_pairs = len(iterator)
        breakpoints = 300
        previous_time = time()

        for counter, (UP1, UP2) in enumerate(iterator):

            idx1, idx2 = (self.inflated_lbl2idx[UP1], self.inflated_lbl2idx[UP2])
            pre_reach = self.UP2GO_Reachable_nodes[UP1] + \
                self.UP2GO_Reachable_nodes[UP2] + [UP1] + [UP2]
            reach = [self.inflated_lbl2idx[label] for label in pre_reach]

            current_upper, voltage_diff = cr.group_edge_current_with_limitations(
                inflated_laplacian=self.inflated_Laplacian,
                idx_pair=(idx1, idx2),
                reach_limiter=reach)

            self.current_accumulator = self.current_accumulator +\
                cr.sparse_abs(current_upper)

            self.UP2UP_voltages[(UP1, UP2)] = voltage_diff

            if counter % breakpoints == 0 and counter > 1:
                # TODO: [load bar] the internal loop load bar goes here
                compops = float(breakpoints) / (time() - previous_time)
                mins_before_termination = (total_pairs-counter) / compops // 60
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
            ln = len(self.active_up_sample)
            self.current_accumulator /= (ln * (ln - 1) / 2)

        if memoized:
            self._dump_memoized()

        index_current = cr.get_current_through_nodes(self.current_accumulator)
        self.node_current = dict((self.inflated_idx2lbl[idx], val)
                                 for idx, val in enumerate(index_current))

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

        for go_term in self.GO2Num.keys():

            if node_current[go_term] > limiting_current:
                characterization_dict[go_term] = [
                    node_current[go_term],
                    self.GO2_Pure_Inf[go_term],
                    len(self.GO2UP_Reachable_nodes[go_term])]

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
            'p_p-value']

        node_char_types = [
            'DOUBLE',
            'VARCHAR',
            'VARCHAR',
            'VARCHAR',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE',
            'DOUBLE']

        if p_value_dict is None:
            p_value_dict = defaultdict(lambda: (np.nan, np.nan, np.nan))

        nan_neg_log10 = lambda x: x if type(x) is str else -np.log10(x)

        char_dict = {}

        if self.sparsely_sampled:
            log.warning('Links between the elements should not be trusted: the computations was '
                        'sampling and was not complete')

        for GO in self.GO2Num.keys():
            char_dict[GO] = [str(self.node_current[GO]),
                             'GO', self.GO_legacy_ids[GO],
                             self.GO_names[GO].replace(',', '-'),
                             str(self.GO2_Pure_Inf[GO]),
                             str(len(self.GO2UP_Reachable_nodes[GO])),
                             str(p_value_dict[int(GO)][0]),
                             str(nan_neg_log10(p_value_dict[int(GO)][0]))]

        for UP in self.active_up_sample:
            char_dict[UP] = [str(self.node_current[UP]),
                             'UP', self.UP_names[UP][0],
                             str(self.UP_names[UP][1]).replace(',', '-'),
                             str(self.binding_intensity),
                             '1',
                             '0.1',
                             '1']  #DOC: document the defaults

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
        gdf_exporter.write()

    def randomly_sample(
            self,
            samples_size,
            samples_each_size,
            sparse_rounds=False,
            memoized=False,
            no_add=False,
            pool_no=None):
        """
        Randomly samples the set of deprecated_reached_uniprots_neo4j_id_list used to create the model.
        This is the null model creation routine

        :param samples_size: list of numbers of uniprots we would like to create the model for
        :param samples_each_size: how many times we would like to sample each unirot number
        :param sparse_rounds:  if we want to use sparse sampling
        (usefull in case of large uniprot sets),
        we would use this option
        :param memoized: if set to True, the sampling would be rememberd for export.
        Usefull in case of the chromosome comparison
        :param no_add: if set to True, the result of sampling will not be added to the database
        of samples. Usefull if re-running tests with similar parameters several times.
        :param pool_no: explicit sampling pool number (used for reporting/debugging)
        :raise Exception: if the number of items in the samples size ann saples_each size are
         different
        """
        if not len(samples_size) == len(samples_each_size):
            raise Exception('Not the same list sizes!')

        # if chromosome_specific:
        #     self.to_deprecate_connected_uniprots = list(set(self.to_deprecate_connected_uniprots).intersection(
        #         set(self.background_set.deprecated_chromosomes_2_uniprot[str(
        #             chromosome_specific)])))

        for sample_size, iterations in zip(samples_size, samples_each_size):

            sample_size = min(sample_size, len(self._background))
            for i in range(0, iterations):
                log.info("selecting %d from %d ups" % (sample_size, len(self._background)))
                shuffle(self._background)
                analytics_up_list = self._background[:sample_size]

                self.set_uniprot_source(analytics_up_list)

                log.info('Sampling thread: %s, Thread hex: %s; '
                         'sampling characteristics: sys_hash: %s, size: %s, '
                         'sparse_rounds: %s' % (pool_no, self.thread_hex,
                                                self.md5_hash(), sample_size,
                                                sparse_rounds))

                self.compute_current_and_potentials(
                    memoized=memoized, sourced=False, sparse_samples=sparse_rounds)

                md5 = hashlib.md5(
                    json.dumps(
                        sorted(analytics_up_list),
                        sort_keys=True).encode('utf-8')).hexdigest()

                if not no_add:
                    log.info("Sampling thread %s: Adding a blanc:"
                             "\t size: %s \t sys_hash: %s \t sparse_rounds: %s, matrix weight: %s" % (
                                pool_no, sample_size, md5, sparse_rounds, np.sum(self.current_accumulator)))

                    insert_annotome_rand_samp(
                        {
                            'UP_hash': md5,
                            'sys_hash': self.md5_hash(),
                            'size': sample_size,
                            'sparse_rounds': sparse_rounds,
                            'UPs': pickle.dumps(analytics_up_list),
                            'currents': pickle.dumps(
                                (self.current_accumulator,
                                 self.node_current)),
                            'voltages': pickle.dumps(
                                self.UP2UP_voltages)})

                if not sparse_rounds:
                    log.info('Sampling thread %s: Thread hex: %s \t'
                             ' Sample size: %s \t iteration: %s\t'
                             ' compop/s: %s \t '
                             'time: %s ',
                             pool_no, self.thread_hex,
                             sample_size, i,
                             "{0:.2f}".format(sample_size * (sample_size - 1) / 2 / self._time()),
                             self.pretty_time())

                else:
                    log.info('Sampling thread %s: Thread hex: %s \t'
                             ' Sample size: %d \t iteration: %d \t'
                             ' compop/s: %s \t '
                             'time: %s, sparse @ %s ',
                             pool_no, self.thread_hex,
                             sample_size, i,
                             "{0:.2f}".format(sample_size * sparse_rounds / 2 / self._time()),
                             self.pretty_time(), sparse_rounds)

                # TODO: [load bar]: the external loop load bar goes here

    def get_independent_linear_groups(self):
        """
        Recovers independent linear groups of the GO terms. Independent linear groups are
        those that share a significant amount of deprecated_reached_uniprots_neo4j_id_list in common
        """
        self.indep_lapl = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
        for GO_list in self.UP2GO_Reachable_nodes.values():
            for GO1, GO2 in combinations(GO_list, 2):
                idx1, idx2 = (self.GO2Num[GO1], self.GO2Num[GO2])
                self.indep_lapl[idx1, idx2] += -1
                self.indep_lapl[idx2, idx1] += -1
                self.indep_lapl[idx2, idx2] += 1
                self.indep_lapl[idx1, idx1] += 1


if __name__ == '__main__':
    # Creates an instance of MatrixGetter and loads pre-computed values

    go_interface_instance = GeneOntologyInterface(
        background_up_ids=get_background_bulbs_ids())
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
    # go_interface_instance.compute_current_and_potentials(sparse_samples=10)
    # go_interface_instance.export_conduction_system()

    # go_interface_instance.deprecated_export_subsystem(experimental, ['186958', '142401', '147798', '164077'])

    # data_array = np.array([log(val) for val in go