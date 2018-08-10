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
from collections import defaultdict
from copy import copy
from csv import reader
from itertools import combinations, chain
from pprint import PrettyPrinter
from random import shuffle
from time import time
import traceback, sys

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import shortest_path

from bioflow.algorithms_bank import conduction_routines as cr
from bioflow.main_configs import Dumps, Outputs, annotome_rand_samp
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.neo4j_db.db_io_routines import get_db_id
from bioflow.utils.gdfExportInterface import GdfExportInterface
from bioflow.utils.io_routines import dump_object, undump_object, get_background_bulbs_ids
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)


def _characterise(_object):
    print 'Object of size %s and type %s' % (len(_object), type(_object))


def _characterise_mat(matrix):
    print 'Matrix of shape %s, type %s and has %s non-zero terms, min is %s, max is %s' % \
          (matrix.shape, type(matrix), len(matrix.nonzero()[0]), '<>', '<>')


class GeneOntologyInterface(object):
    """
    General class to recover all the information associated with GO from database and buffer
     them for further use.

    :param Filter:
    :param Uniprot_Node_IDs: A list of hte reached_uniprots_neo4j_id_list that will be used
     for the GO reach and informativity computation. Beyond database and formalism issues,
     this allows to adapt method when limited list of UP is of interest
    :param correction_factor:
    :param Ultraspec_clean:
    :param Ultraspec_lvl: parameter how much uniprots have to point to a GO term for it not
     to be considered ultra-specific anymore
    """

    _GOUpTypes = ["is_a_go", "is_part_of_go"]
    _GORegTypes = ["is_Regulant"]

    def __init__(self, namespace_filter=('biological_process',),
                 uniprot_node_ids=None,
                 correction_factor=(1, 1),
                 ultraspec_clean=True, ultraspec_lvl=3):

        self.interactome_interface_instance = InteractomeInterface(True, True)
        self.interactome_interface_instance.fast_load()
        init_set = self.interactome_interface_instance.all_uniprots_neo4j_id_list

        if uniprot_node_ids:
            init_set = list(set(init_set).intersection(set(uniprot_node_ids)))

        self.go_namespace_filter = list(namespace_filter)
        self.InitSet = init_set
        self.correction_factor = correction_factor
        self.ultraspec_cleaned = ultraspec_clean
        self.ultraspec_lvl = ultraspec_lvl
        self.init_time = time()
        self.partial_time = time()

        self.UPs_without_GO = set()
        self.UP2GO_Dict = {}
        self.GO2UP = defaultdict(list)
        self.SeedSet = set()
        self.All_GOs = []
        self.GO2Num = {}
        self.Num2GO = {}
        self.total_Entropy = None

        self.Reachable_nodes_dict = {}

        self.UP2GO_Reachable_nodes = {}
        self.GO2UP_Reachable_nodes = {}
        self.UP2GO_step_Reachable_nodes = {}
        self.GO2UP_step_Reachable_nodes = {}
        self.GO2_Pure_Inf = {}
        self.GO2_Weighted_Ent = {}

        self.GO_Names = {}
        self.GO_Legacy_IDs = {}
        self.rev_GO_IDs = {}

        self.UP_Names = {}

        self.adjacency_matrix = np.zeros((2, 2))
        self.dir_adj_matrix = np.zeros((2, 2))
        self.laplacian_matrix = np.zeros((2, 2))
        self.weighted_laplacian_matrix = np.zeros((2, 2))

        self.Sign_retaining_matrix = np.zeros((2, 2))

        self.TimesReached = {}
        self.accelerationDict = {}
        self.Reverse_Dict = {}
        self.GO_Node_ID2Reach = {}

        # everytghing below should not be dumped, but exported to the
        # mongoDB

        self.inflated_Laplacian = np.zeros((2, 2))
        self.inflated_idx2lbl = {}
        self.inflated_lbl2idx = {}
        self.binding_intensity = 0

        self.analytic_uniprots = []
        self.UP2UP_voltages = {}
        self.uniprots_2_voltage_and_circulation = {}  # can be safely renamed

        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}
        self.call_coutner = 0

        char_set = string.ascii_uppercase + string.digits
        self.random_tag = ''.join(random.sample(char_set * 6, 6))

        self.Indep_Lapl = np.zeros((2, 2))
        self.uncomplete_compute = False

        self.main_set = self.InitSet

        log.info('Setting up GO Interface with namespaces %s and %s root UPs',
                 self.go_namespace_filter, len(self.InitSet))

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

    def dump_statics(self):
        dump_object(
            Dumps.GO_builder_stat,
            (self.go_namespace_filter,
             self.InitSet,
             self.correction_factor,
             self.ultraspec_cleaned,
             self.ultraspec_lvl))

    @staticmethod
    def undump_statics():
        return undump_object(Dumps.GO_builder_stat)

    def dump_core(self):
        dump_object(
            Dumps.GO_dump,
            (self.UP2GO_Dict,
             self.GO2UP,
             self.SeedSet,
             self.Reachable_nodes_dict,
             self.GO_Names,
             self.GO_Legacy_IDs,
             self.rev_GO_IDs,
             self.All_GOs,
             self.GO2Num,
             self.Num2GO,
             self.UP_Names,
             self.UPs_without_GO))

    def undump_core(self):
        self.UP2GO_Dict, self.GO2UP, self.SeedSet, self.Reachable_nodes_dict,\
            self.GO_Names, self.GO_Legacy_IDs, self.rev_GO_IDs, self.All_GOs,\
            self.GO2Num, self.Num2GO, self.UP_Names, self.UPs_without_GO =\
            undump_object(Dumps.GO_dump)

    def dump_matrices(self):
        dump_object(
            Dumps.GO_Mats,
            (self.adjacency_matrix,
             self.dir_adj_matrix,
             self.laplacian_matrix))

    def undump_matrices(self):
        self.adjacency_matrix, self.dir_adj_matrix, self.laplacian_matrix = undump_object(
            Dumps.GO_Mats)

    def dump_informativities(self):
        dump_object(
            Dumps.GO_Infos,
            (self.UP2GO_Reachable_nodes,
             self.GO2UP_Reachable_nodes,
             self.UP2GO_step_Reachable_nodes,
             self.GO2UP_step_Reachable_nodes,
             self.GO2_Pure_Inf,
             self.GO2_Weighted_Ent))

    def undump_informativities(self):
        self.UP2GO_Reachable_nodes, self.GO2UP_Reachable_nodes, self.UP2GO_step_Reachable_nodes, \
            self.GO2UP_step_Reachable_nodes, self.GO2_Pure_Inf, self.GO2_Weighted_Ent = \
            undump_object(Dumps.GO_Infos)

    def dump_inflated_elements(self):
        dump_object(
            Dumps.GO_Inflated,
            (self.inflated_Laplacian,
             self.inflated_idx2lbl,
             self.inflated_lbl2idx,
             self.binding_intensity))

    def undump_inflated_elements(self):
        self.inflated_Laplacian, self.inflated_idx2lbl, \
            self.inflated_lbl2idx, self.binding_intensity = \
            undump_object(Dumps.GO_Inflated)

    def dump_memoized(self):
        md5 = hashlib.md5(
            json.dumps(
                sorted(
                    self.analytic_uniprots),
                sort_keys=True)).hexdigest()
        payload = {
            'UP_hash': md5, 'sys_hash': self.md5_hash(), 'size': len(
                self.analytic_uniprots), 'UPs': pickle.dumps(
                self.analytic_uniprots), 'currents': pickle.dumps(
                (self.current_accumulator, self.node_current)), 'voltages': pickle.dumps(
                    self.uniprots_2_voltage_and_circulation)}
        dump_object(Dumps.GO_Analysis_memoized, payload)

    @staticmethod
    def undump_memoized():
        """
        :return: undumped memoized analysis
        """
        return undump_object(Dumps.GO_Analysis_memoized)

    def dump_independent_linear_sets(self):
        dump_object(Dumps.GO_Indep_Linset, self.Indep_Lapl)

    def undump_independent_linear_sets(self):
        self.Indep_Lapl = undump_object(Dumps.GO_Indep_Linset)

    def full_rebuild(self):
        self.get_gene_ontology_access()
        self.get_gene_ontology_structure()
        self.get_go_adjacency_and_laplacian()
        self.get_go_reach()
        if self.ultraspec_cleaned:
            self.filter_out_too_specific()
        self.get_laplacians()
        self.inflate_matrix_and_indexes()

        self.dump_statics()
        self.dump_core()
        self.dump_matrices()
        self.dump_informativities()
        self.dump_inflated_elements()

        log.info('Finished rebuilding the GO Interface object %s', self.pretty_time())

    def load(self):
        """
        loads itself from the saved dumps, in case the Filtering system is the same

        """
        namespace_filter, initial_set, correction_factor, ultraspec_cleaned, ultraspec_lvl = \
            self.undump_statics()
        if self.go_namespace_filter != namespace_filter:
            log.critical("Wrong Filtering attempted to be recovered from storage")
            raise Exception(
                "Wrong Filtering attempted to be recovered from storage")
        if self.InitSet != initial_set:
            print len(self.InitSet)
            print len(initial_set)
            print traceback.print_stack()
            log.critical("Wrong initial_set attempted to be recovered from storage")
            raise Exception(
                "Wrong initial_set attempted to be recovered from storage")
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
        self.undump_core()
        self.undump_matrices()
        self.undump_informativities()
        self.undump_inflated_elements()

    def get_gene_ontology_access(self):
        """
        Loads all of the relations between the UNIPROTs and GOs as one giant dictionary

        """
        uniprots_without_gene_ontology_terms = 0
        log.info('Starting GO matrix mapping starting from %s uniprots', len(self.InitSet))
        for uniprot_neo4j_id in self.InitSet:
            uniprot_specific_gos = []

            up_node = DatabaseGraph.get(uniprot_neo4j_id)
            self.UP_Names[uniprot_neo4j_id] = [up_node.properties['legacyId'],
                                               up_node.properties['displayName']]
            attached_go_nodes = DatabaseGraph.get_linked(uniprot_neo4j_id,
                                                         link_type='is_go_annotation')
            for go_node in attached_go_nodes:
                if go_node.properties['Namespace'] in self.go_namespace_filter:
                    go_node_neo4j_id = get_db_id(go_node)
                    uniprot_specific_gos.append(go_node_neo4j_id)
                    self.GO2UP[go_node_neo4j_id].append(uniprot_neo4j_id)
                    self.SeedSet.add(go_node_neo4j_id)

            if not uniprot_specific_gos:
                uniprots_without_gene_ontology_terms += 1
                log.debug("UP without GO was found. UP bulbs_id: %s, \t name: %s",
                          uniprot_neo4j_id, self.UP_Names[uniprot_neo4j_id])
                self.UPs_without_GO.add(uniprot_neo4j_id)
            else:
                self.UP2GO_Dict[uniprot_neo4j_id] = copy(uniprot_specific_gos)

        log.info('total number of UPs without a go_node annotation: %s out of %s',
                 uniprots_without_gene_ontology_terms, len(self.InitSet))

    # TODO: REFACTORING. Method is excessively complex.
    def get_gene_ontology_structure(self):
        """
        Loads all of the relations between the GOs that are generalisation of the seedList
         GOs and that are withing the types specified in go_namespace_filter

        """
        visited_set = set()
        seeds_list = copy(list(self.SeedSet))
        log.info('Starting gene ontology structure retrieval from the set of %s seeds',
                 len(self.SeedSet))

        while seeds_list:
            node_id = seeds_list.pop()
            visited_set.add(node_id)
            local_uniprot_list = []
            local_regulation_list = []
            local_up_regulation_list = []
            local_down_regulation_list = []
            gene_ontology_node = DatabaseGraph.get(node_id, 'GOTerm')
            self.GO_Names[node_id] = str(gene_ontology_node.properties['displayName'])
            self.GO_Legacy_IDs[node_id] = str(gene_ontology_node.properties['legacyId'])
            self.rev_GO_IDs[gene_ontology_node.properties['legacyId']] = node_id

            for relation_type in chain(self._GOUpTypes, self._GORegTypes):
                related_go_nodes = DatabaseGraph.get_linked(node_id, 'out', relation_type)

                if not related_go_nodes:
                    continue  # skip in case GO Node has no outgoing relations to other GO nodes
                for go_node in related_go_nodes:
                    if go_node.properties['Namespace'] not in self.go_namespace_filter:
                        continue
                    node_bulbs_id = get_db_id(go_node)
                    if node_bulbs_id not in visited_set:
                        seeds_list.append(node_bulbs_id)
                    if relation_type in self._GOUpTypes:
                        local_uniprot_list.append(node_bulbs_id)
                    else:
                        local_regulation_list.append(node_bulbs_id)

                rev_generator = DatabaseGraph.get_linked(node_id, 'in', relation_type)

                if not rev_generator:
                    continue
                for go_node in rev_generator:
                    if go_node.properties['Namespace'] not in self.go_namespace_filter:
                        continue
                    node_bulbs_id = get_db_id(go_node)
                    if relation_type in self._GOUpTypes:
                        local_down_regulation_list.append(node_bulbs_id)
                    else:
                        local_up_regulation_list.append(node_bulbs_id)

            self.Reachable_nodes_dict[node_id] = (
                list(set(local_uniprot_list)),
                list(set(local_regulation_list)),
                list(set(local_down_regulation_list)),
                list(set(local_up_regulation_list)))

        self.All_GOs = list(visited_set)
        self.Num2GO = dict((i, val) for i, val in enumerate(self.All_GOs))
        self.GO2Num = dict((val, i) for i, val in enumerate(self.All_GOs))

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
            for node, package in self.Reachable_nodes_dict.iteritems():
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
            for node, package in self.Reachable_nodes_dict.iteritems():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.GO2Num[node], self.GO2Num[node2])
                    base_matrix[idx] = 1

            self.dir_adj_matrix = copy(base_matrix)

        build_adjacency()
        build_dir_adj()

    def calculate_informativity(self, number):
        """
        returns an entropy given by a number of equi-probable events, where event is the number.

        :param number:
        """
        if number < 1.0:
            log.critical("Wrong value provided for entropy computation")
            raise Exception("Wrong value provided for entropy computation")
        if not self.total_Entropy:
            self.total_Entropy = - \
                math.log(1 / float(len(self.UP2GO_Dict.keys())), 2)
        if number == 1.0:
            return 2 * self.total_Entropy
        return pow(-self.correction_factor[0] * self.total_Entropy /
                   math.log(1 / float(number), 2), self.correction_factor[1])

    # TODO: REFACTORING. Method is excessively complex.
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
            dict_len = {key: [len(val), len(step_reach[key].keys())]
                        for key, val in reach.iteritems()}
            for key, val in dict_len.iteritems():
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
            for key, val_list in _val_dict.iteritems():
                summer += filter_function(key) * len(val_list)
            return summer

        dir_reg_path = shortest_path(
            self.dir_adj_matrix, directed=True, method='D')
        dir_reg_path[np.isinf(dir_reg_path)] = 0.0
        dir_reg_path = lil_matrix(dir_reg_path)

        self.GO2UP_Reachable_nodes = dict(
            (el, []) for el in self.Reachable_nodes_dict.keys())
        self.GO2UP_Reachable_nodes.update(self.GO2UP)

        pre_go2up_step_reachable_nodes = \
            dict((key, dict((v, 0) for v in val))
                 for key, val in self.GO2UP_Reachable_nodes.iteritems())
        # when called on possibly un-encoutenred items, anticipate a default
        # value of 10 000

        # Now just scan vertical columns and add UP terms attached
        for idx1, idx2 in zip(list(dir_reg_path.nonzero()[0]),
                              list(dir_reg_path.nonzero()[1])):
            self.GO2UP_Reachable_nodes[self.Num2GO[idx2]] +=\
                self.GO2UP_Reachable_nodes[self.Num2GO[idx1]]
            if dir_reg_path[idx1, idx2] < 1.0:
                log.critical("null in non-null patch")
                raise Exception("null in non-null patch")
            step_reach_upgrade = dict(
                (key,
                 val +
                 dir_reg_path[
                     idx1,
                     idx2]) for key,
                val in pre_go2up_step_reachable_nodes[
                    self.Num2GO[idx1]].iteritems())
            for k, v in step_reach_upgrade.iteritems():
                pre_go2up_step_reachable_nodes[
                    self.Num2GO[idx2]][k] = min(
                    pre_go2up_step_reachable_nodes[
                        self.Num2GO[idx2]].setdefault(
                        k, 100000), v)

        for key, val in self.GO2UP_Reachable_nodes.iteritems():
            self.GO2UP_Reachable_nodes[key] = list(set(val))

        verify_equivalence_of_reaches(
            pre_go2up_step_reachable_nodes,
            self.GO2UP_Reachable_nodes)

        # Now we need to invert the reach to get the set of all the primary and
        # derived GO terms that describe a UP
        self.UP2GO_Reachable_nodes = dict(
            (key, []) for key in self.UP2GO_Dict.keys())
        self.UP2GO_step_Reachable_nodes = dict(
            (key, defaultdict(list)) for key in self.UP2GO_Dict.keys())
        self.GO2UP_step_Reachable_nodes = dict(
            (key, defaultdict(list)) for key in pre_go2up_step_reachable_nodes.keys())
        for key, val_dict in pre_go2up_step_reachable_nodes.iteritems():
            for k, v in val_dict.iteritems():
                self.GO2UP_step_Reachable_nodes[key][v].append(k)
                self.UP2GO_step_Reachable_nodes[k][v].append(key)
                self.UP2GO_Reachable_nodes[k].append(key)

        # and finally we compute the pure and weighted informativity for each
        # term
        self.GO2_Pure_Inf = dict((key, self.calculate_informativity(len(val)))
                                 for key, val in self.GO2UP_Reachable_nodes.iteritems())
        self.GO2_Weighted_Ent = dict(
            (key,
             self.calculate_informativity(special_sum(val_dict)))
            for key, val_dict in self.GO2UP_step_Reachable_nodes.iteritems())

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
            zip(list(base_matrix.nonzero()[0]), list(base_matrix.nonzero()[1])))

        for idx1, idx2 in nz_list:
            min_inf = min(
                self.GO2_Pure_Inf[
                    self.Num2GO[idx1]], self.GO2_Pure_Inf[
                    self.Num2GO[idx2]])
            base_matrix[idx1, idx2] = -min_inf
            base_matrix[idx2, idx1] = -min_inf
            base_matrix[idx2, idx2] += min_inf
            base_matrix[idx1, idx1] += min_inf

        self.laplacian_matrix = base_matrix

    def compute_uniprot_dict(self):
        """
        Computes the uniprot method required by some other dictionary

        :return:
        """
        uniprot_dict = {}
        for elt in self.interactome_interface_instance.reached_uniprots_neo4j_id_list:
            node = DatabaseGraph.get(elt, 'UNIPROT')
            alt_id = node.properties['legacyId']
            # TODO: now can be suppressed
            uniprot_dict[alt_id] = (
                elt, self.interactome_interface_instance.neo4j_id_2_display_name[elt])
            uniprot_dict[elt] = alt_id
        pickle.dump(uniprot_dict, file(Dumps.Up_dict_dump, 'w'))
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
                                  in self.GO2UP_Reachable_nodes.iteritems()
                                  if len(reach) < self.ultraspec_lvl)
        for GO in ultraspec_go_terms:
            self.GO2_Pure_Inf[GO] = rep_val

    def md5_hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been
        correctly defined before dump/retrieval
        """
        sorted_initset = sorted(self.InitSet)
        data = [
            self.go_namespace_filter,
            sorted_initset,
            self.correction_factor,
            self.ultraspec_cleaned,
            self.ultraspec_lvl]
        md5 = hashlib.md5(json.dumps(data, sort_keys=True)).hexdigest()
        return str(md5)

    def inflate_matrix_and_indexes(self):
        """
        Performs the laplacian matrix inflation to incorporate the uniprots on which we
         will be running the
        """
        # branching distribution: at least 10x the biggest conductivity of the
        # system, unless too specific, in which case ~ specific level
        self.binding_intensity = 10 * self.calculate_informativity(self.ultraspec_lvl)
        fixed_index = self.laplacian_matrix.shape[0]
        self_connectable_uniprots = list(
            set(self.InitSet) - set(self.UPs_without_GO))
        up2idxs = dict((UP, fixed_index + Idx)
                       for Idx, UP in enumerate(self_connectable_uniprots))
        idx2ups = dict((Idx, UP) for UP, Idx in up2idxs.iteritems())
        self.inflated_Laplacian = lil_matrix(
            (self.laplacian_matrix.shape[0] + len(self_connectable_uniprots),
             self.laplacian_matrix.shape[1] + len(self_connectable_uniprots)))
        self.inflated_Laplacian[:self.laplacian_matrix.shape[0], :self.laplacian_matrix.shape[1]] =\
            self.laplacian_matrix

        for uniprot in self_connectable_uniprots:
            for go_term in self.UP2GO_Dict[uniprot]:
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
        Sets the reached_uniprots_neo4j_id_list on which the circulation computation routines
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
        self.analytic_uniprots = [
            uniprot for uniprot in uniprots if uniprot in self.UP2GO_Dict.keys()]

    def build_extended_conduction_system(
            self,
            memoized=True,
            sourced=False,
            incremental=False,
            cancellation=True,
            sparse_samples=False):
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier
        knowledge flow analysis


        :param memoized: if the tensions and individual relation matrices should be stored in
         the matrix and dumped at the end computation (required for submatrix re-computation)
        :param sourced: if true, all the relations will be looked up and not computed. Useful
        for the retrieval of sub-circulation group, but requires the
        uniprots_2_voltage_and_circulation to be pre-filled
        :param incremental: if True, all the circulation computation will be added to the
        existing ones. Useful for the computation of particularly big systems with
        intermediate dumps
        :param cancellation: divides the final current by #Nodes**2/2, i.e. makes the currents
        comparable between circulation systems of different  sizes.
        :param sparse_samples: if set to an integer the sampling will be sparse and not dense,
         i.e. instead of computation for each node pair, only an estimation will be made, equal to
         computing sparse_samples association with other randomly chosen nodes
        :type sparse_samples: int
        :return: adjusted conduction system
        """
        if not incremental or self.current_accumulator == np.zeros((2, 2)):
            self.current_accumulator = lil_matrix(self.inflated_Laplacian.shape)
            self.UP2UP_voltages = {}
            if not sourced:
                self.uniprots_2_voltage_and_circulation = {}

        iterator = []
        if sparse_samples:
            for _ in range(0, sparse_samples):
                _length = copy(self.analytic_uniprots)
                random.shuffle(_length)
                iterator += zip(_length[:len(_length) / 2], _length[len(_length) / 2:])
                self.uncomplete_compute = True
        else:
            iterator = combinations(self.analytic_uniprots, 2)

        iterator = [item for item in iterator]

        total_pairs = len(iterator)
        breakpoints = 300
        previous_time = time()

        for counter, (UP1, UP2) in enumerate(iterator):

            if sourced:
                self.current_accumulator = self.current_accumulator + \
                    cr.sparse_abs(self.uniprots_2_voltage_and_circulation[
                                  tuple(sorted((UP1, UP2)))][1])
                continue

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

            if memoized:
                self.uniprots_2_voltage_and_circulation[
                    tuple(sorted((UP1, UP2)))] = \
                    (voltage_diff, current_upper)

            if counter % breakpoints == 0 and counter > 1:
                compops = float(breakpoints) / (time() - previous_time)
                log.info("progress: %s/%s, current speed: %s compops, time remaining: %s min"
                         % (counter, total_pairs, compops, (total_pairs - counter) / compops / 60))
                previous_time = time()

        if cancellation:  # TODO: factor that one into the Conduction Routilens
            ln = len(self.analytic_uniprots)
            self.current_accumulator /= (ln * (ln - 1) / 2)

        if memoized:
            self.dump_memoized()

        index_current = cr.get_current_through_nodes(self.current_accumulator)
        self.node_current = dict(
            (self.inflated_idx2lbl[idx],
             val) for idx,
            val in enumerate(index_current))

    def format_node_props(self, node_current, limit=0.01):
        """
        Formats the nodes for the analysis by in the knowledge_access_analysis module

        :param node_current: Current through the GO nodes
        :param limit: hard limit to go_namespace_filter out the GO terms with too little current
         (compensates the minor currents in the gird)
        :return: {GO:[node current, pure GO informativity, Number of reachable nodes]}
        """
        char_dict = {}
        limiting_current = max(node_current.values()) * limit
        for go_term in self.GO2Num.iterkeys():

            if node_current[go_term] > limiting_current:
                char_dict[go_term] = [
                    node_current[go_term],
                    self.GO2_Pure_Inf[go_term],
                    len(self.GO2UP_Reachable_nodes[go_term])]

        return char_dict

    def export_conduction_system(self):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format
         and flushes it into a file that can be viewed with Gephi

        :raise Warning:
        """
        node_char_names = [
            'Current',
            'Type',
            'Legacy_ID',
            'Names',
            'Pure_informativity',
            'Confusion_potential']
        node_char_types = [
            'DOUBLE',
            'VARCHAR',
            'VARCHAR',
            'VARCHAR',
            'DOUBLE',
            'DOUBLE']
        char_dict = {}

        if self.uncomplete_compute:
            log.warning('Links between the elements should not be trusted: the computations was '
                        'sampling and was not complete')

        for GO in self.GO2Num.iterkeys():
            char_dict[GO] = [str(self.node_current[GO]),
                             'GO', self.GO_Legacy_IDs[GO],
                             self.GO_Names[GO].replace(',', '-'),
                             str(self.GO2_Pure_Inf[GO]),
                             str(len(self.GO2UP_Reachable_nodes[GO]))]

        for UP in self.analytic_uniprots:
            char_dict[UP] = [str(self.node_current[UP]),
                             'UP', self.UP_Names[UP][0],
                             str(self.UP_Names[UP][1]).replace(',', '-'),
                             str(self.binding_intensity),
                             '1']

        gdf_exporter = GdfExportInterface(
            target_fname=Outputs.GO_GDF_output,
            field_names=node_char_names,
            field_types=node_char_types,
            node_properties_dict=char_dict,
            min_current=0.01,
            index_2_label=self.inflated_idx2lbl,
            label_2_index=self.inflated_lbl2idx,
            current_matrix=self.current_accumulator)
        gdf_exporter.write()

    def export_subsystem(self, uniprot_system, uniprot_subsystem):
        """
        Exports the subsystem of reached_uniprots_neo4j_id_list and circulation between
        them based on a larger precalculated system.This is possible only of the memoization
        parameter was on during the execution of "build_extended_circulation_system()"
        function execution.

        :param uniprot_system: The set of uniprots for which the larger system was calculated
        :param uniprot_subsystem: the set of reached_uniprots_neo4j_id_list we are interested in
        :raise Exception: if the set of uniprots for which the larger system was calculated
         doesn't correspond to what is stored in the dumps
        """
        current_recombinator = self.undump_memoized()
        if not set(uniprot_system) == set(
                pickle.loads(current_recombinator['UPs'])):
            raise Exception('Wrong UP system re-analyzed')
        self.uniprots_2_voltage_and_circulation = pickle.loads(
            current_recombinator['voltages'])
        self.set_uniprot_source(uniprot_subsystem)
        self.build_extended_conduction_system(memoized=False, sourced=True)
        self.export_conduction_system()

    def randomly_sample(
            self,
            samples_size,
            samples_each_size,
            sparse_rounds=False,
            chromosome_specific=False,
            memoized=False,
            no_add=False):
        """
        Randomly samples the set of reached_uniprots_neo4j_id_list used to create the model.

         This is the null model creation routine

        :param samples_size: list of numbers of uniprots we would like to create the model for
        :param samples_each_size: how many times we would like to sample each unirot number
        :param sparse_rounds:  if we want to use sparse sampling
        (usefull in case of large uniprot sets),
        we would use this option
        :param chromosome_specific: if we want the sampling to be chromosome-specific,
        set this parameter to the
        number of chromosome to sample from
        :param memoized: if set to True, the sampling would be rememberd for export.
        Usefull in case of the chromosome comparison
        :param no_add: if set to True, the result of sampling will not be added to the database
        of samples. Usefull if re-running tests with similar parameters several times.
        :raise Exception: if the number of items in the samples size ann saples_each size are
         different
        """
        if not len(samples_size) == len(samples_each_size):
            raise Exception('Not the same list sizes!')

        self_connectable_uniprots = list(
            set(self.InitSet) - set(self.UPs_without_GO))

        if chromosome_specific:
            self_connectable_uniprots = list(set(self_connectable_uniprots).intersection(
                set(self.interactome_interface_instance.chromosomes_2_uniprot[str(
                    chromosome_specific)])))

        for sample_size, iterations in zip(samples_size, samples_each_size):
            sample_size = min(sample_size, len(self_connectable_uniprots))
            for i in range(0, iterations):
                shuffle(self_connectable_uniprots)
                analytics_up_list = self_connectable_uniprots[:sample_size]
                self.set_uniprot_source(analytics_up_list)
                self.build_extended_conduction_system(
                    memoized=memoized, sourced=False, sparse_samples=sparse_rounds)

                md5 = hashlib.md5(
                    json.dumps(
                        sorted(analytics_up_list),
                        sort_keys=True)).hexdigest()

                if not no_add:
                    annotome_rand_samp.insert(
                        {
                            'UP_hash': md5,
                            'sys_hash': self.md5_hash(),
                            'size': sample_size,
                            'chrom': str(chromosome_specific),
                            'sparse_rounds': sparse_rounds,
                            'UPs': pickle.dumps(analytics_up_list),
                            'currents': pickle.dumps(
                                (self.current_accumulator,
                                 self.node_current)),
                            'voltages': pickle.dumps(
                                self.UP2UP_voltages)})

                if not sparse_rounds:

                    log.info(
                        'Random ID: %s \t Sample size: %s \t iteration: %s\t compops: %s \t time: %s ',
                        self.random_tag, sample_size, i,
                        "{0:.2f}".format(sample_size * (sample_size - 1) / 2 / self._time()),
                        self.pretty_time())

                else:
                    log.info('Random ID: %s \t Sample size: %s \t iteration: %s\t compops: %s \t '
                             'time: %s ', self.random_tag, sample_size, i,
                             "{0:.2f}".format(sample_size * sparse_rounds / 2 / self._time()),
                             self.pretty_time())

    def get_independent_linear_groups(self):
        """
        Recovers independent linear groups of the GO terms. Independent linear groups are
        those that share a significant amount of reached_uniprots_neo4j_id_list in common
        """
        self.Indep_Lapl = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
        for GO_list in self.UP2GO_Reachable_nodes.itervalues():
            for GO1, GO2 in combinations(GO_list, 2):
                idx1, idx2 = (self.GO2Num[GO1], self.GO2Num[GO2])
                self.Indep_Lapl[idx1, idx2] += -1
                self.Indep_Lapl[idx2, idx1] += -1
                self.Indep_Lapl[idx2, idx2] += 1
                self.Indep_Lapl[idx1, idx1] += 1


if __name__ == '__main__':
    # Creates an instance of MatrixGetter and loads pre-computed values

    go_interface_instance = GeneOntologyInterface(uniprot_node_ids=get_background_bulbs_ids())
    go_interface_instance.full_rebuild()

    # loading takes 1-6 seconds.
    # fill for reach only is done in 2 seconds,
    # tepping takes another 15,
    # inverting + info computation - 1 more second
    # Laplacian building =>
    ##
    # full computation - 3 minutes 18 seconds; save 7 seconds, retrieval - 3
    # seconds

    # go_interface_instance.load()
    # print go_interface_instance.pretty_time()

    # go_interface_instance.get_indep_linear_groups()
    # go_interface_instance.dump_Indep_Linset()

    # go_interface_instance.randomly_sample([10, 25], [5]*2, chromosome_specific=15)

    # go_interface_instance.set_Uniprot_source(experimental)
    # go_interface_instance.build_extended_conduction_system(sparse_samples=10)
    # go_interface_instance.export_conduction_system()

    # go_interface_instance.export_subsystem(experimental, ['186958', '142401', '147798', '164077'])

    # data_array = np.array([log(val) for val in go_interface_instance.GO2_Pure_Inf.itervalues()])
    # hist(data_array, 100, log=True)
    # show()
