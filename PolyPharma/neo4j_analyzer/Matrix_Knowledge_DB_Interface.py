__author__ = 'ank'
"""
Contains all the tools necessary to map GO ontology and Pathway classification from the database to an Adjacency and
Laplacian graph.
"""

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
import copy
from time import time
import pickle
import operator
from random import shuffle
import math
from scipy.stats.kde import  gaussian_kde
import numpy as np
from pylab import plot, hist, show
from itertools import combinations
from pprint import PrettyPrinter
from math import log
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import shortest_path, connected_components
from scipy.sparse.csgraph._validation import validate_graph
from collections import defaultdict
from itertools import chain
from PolyPharma.neo4j_analyzer.Matrix_Interactome_DB_interface import MatrixGetter
from PolyPharma.configs import Dumps, Outputs, UP_rand_samp
from PolyPharma.neo4j_analyzer.knowledge_access import acceleratedInsert
from PolyPharma.neo4j_analyzer.IO_Routines import dump_object, write_to_csv, undump_object
from PolyPharma.Utils.GDF_export import GDF_export_Interface
import PolyPharma.neo4j_analyzer.Conduction_routines as CR
import hashlib
import json
import random
import string
from multiprocessing import Pool


# Creates an instance of MatrixGetter and loads pre-computed values
MG = MatrixGetter(True, False)
MG.fast_load()


# TODO: switch to the usage of Uniprot set that is independent from the Matrix_Getter, but instead is supplide by the user
    # MG.Uniprot is just an option, even though a very importatn one

# specify the relations that lead to a more general or to an equally regulated node.
GOUpTypes = ["is_a_go", "is_part_of_go"]
GORegTypes = ["is_Regulant"]


ppritner = PrettyPrinter(indent = 4)


def _characterise(objekt):
    print 'Object of size %s and type %s' %(len(objekt),type(objekt))


def _characterise_mat(matrix):
    print 'Matrix of shape %s, type %s and has %s non-zero terms, min is %s, max is %s' %(matrix.shape, type(matrix),
                                                                len(matrix.nonzero()[0]), '<>', '<>')


class GO_Interface(object):
    """
    General calss to recover all the informations associated with GO from database and buffer them for further use.

    :param Filter:
    :param Uniprot_Node_IDs: A list of hte Uniprots that will be used for the GO reach and informativity computation. Beyond
                                database and formalism issues, this allows to adapt method when limited list of UP is of interest
    :param corrfactor:
    :param Ultraspec_clean:
    :param Ultraspec_lvl: parameter how much uniprots have to point to a GO term for it not to be considered ultraspecific anymore
    """

    def __init__(self, Filter, Uniprot_Node_IDs, corrfactor, Ultraspec_clean = False, Ultraspec_lvl = 3,):
        self.Filtr = Filter
        self.InitSet = Uniprot_Node_IDs
        self.corrfactor = corrfactor
        self.ultraspec_cleaned = Ultraspec_clean
        self.ultraspec_lvl = Ultraspec_lvl
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

        self.Reachable_nodes_dict ={}

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

        self.Adjacency_matrix = np.zeros((2, 2))
        self.dir_adj_matrix = np.zeros((2, 2))
        self.Laplacian_matrix = np.zeros((2, 2))
        self.Weighted_Laplacian_matrix = np.zeros((2, 2))

        self.Sign_retaining_matrix = np.zeros((2, 2))

        self.TimesReached = {}
        self.accelerationDict = {}
        self.Reverse_Dict = {}
        self.GO_Node_ID2Reach = {}

        # everytghing below should not be dumped, but should be exported to the mongoDB

        self.inflated_Laplacian = np.zeros((2,2))
        self.inflated_idx2lbl = {}
        self.inflated_lbl2idx = {}
        self.binding_intesity = 0

        self.analytic_Uniprots = []
        self.UP2UP_voltages ={}
        self.UP2circ_and_voltage = {}

        self.current_accumulator = np.zeros((2,2))
        self.node_current = {}
        self.call_coutner = 0

        char_set = string.ascii_uppercase + string.digits
        self.r_ID = ''.join(random.sample(char_set*6, 6))

        self.Indep_Lapl = np.zeros((2,2))


    def pretty_time(self):
        """
        Times the execution

        :return: tuple containing the time since the creation of the Matrix_getter object and since the last cal of function formatted as string
        :rtype: str
        """
        it, pt = (round(time() - self.init_time), round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (int(it) / 60, it % 60, int(pt) / 60, pt % 60)
        self.partial_time = time()
        return pload


    def _time(self):
        pt = time() - self.partial_time
        return pt


    def call_show(self):
        self.call_coutner +=1
        print '*',
        if self.call_coutner % 100 == 99:
            print ''


    def dump_statics(self):
        dump_object(Dumps.GO_builder_stat, (self.Filtr, self.InitSet ,self.corrfactor, self.ultraspec_cleaned, self.ultraspec_lvl))


    def undump_statics(self):
        return undump_object(Dumps.GO_builder_stat)


    def dump_core(self):
        dump_object(Dumps.GO_dump, (self.UP2GO_Dict, self.GO2UP, self.SeedSet, self.Reachable_nodes_dict, self.GO_Names,
                                    self.GO_Legacy_IDs, self.rev_GO_IDs, self.All_GOs, self.GO2Num, self.Num2GO,
                                    self.UP_Names, self.UPs_without_GO))


    def undump_core(self):
        self.UP2GO_Dict, self.GO2UP, self.SeedSet, self.Reachable_nodes_dict, self.GO_Names, self.GO_Legacy_IDs,\
        self.rev_GO_IDs, self.All_GOs, self.GO2Num, self.Num2GO, self.UP_Names, self.UPs_without_GO = undump_object(Dumps.GO_dump)


    def dump_matrices(self):
        dump_object(Dumps.GO_Mats, (self.Adjacency_matrix, self.dir_adj_matrix, self.Laplacian_matrix))


    def undump_matrices(self):
        self.Adjacency_matrix, self.dir_adj_matrix, self.Laplacian_matrix = undump_object(Dumps.GO_Mats)


    def dump_informativities(self):
        dump_object(Dumps.GO_Infos, (self.UP2GO_Reachable_nodes, self.GO2UP_Reachable_nodes, self.UP2GO_step_Reachable_nodes,
                                    self.GO2UP_step_Reachable_nodes, self.GO2_Pure_Inf, self.GO2_Weighted_Ent))


    def undump_informativities(self):
        self.UP2GO_Reachable_nodes, self.GO2UP_Reachable_nodes, self.UP2GO_step_Reachable_nodes, \
        self.GO2UP_step_Reachable_nodes, self.GO2_Pure_Inf, self.GO2_Weighted_Ent = undump_object(Dumps.GO_Infos)


    def dump_inflated_elements(self):
        dump_object(Dumps.GO_Inflated, (self.inflated_Laplacian, self.inflated_idx2lbl, self.inflated_lbl2idx, self.binding_intesity))


    def undump_inflated_elements(self):
        self.inflated_Laplacian, self.inflated_idx2lbl, self.inflated_lbl2idx, self.binding_intesity = undump_object(Dumps.GO_Inflated)


    def dump_memoized(self):
        md5 = hashlib.md5(json.dumps( sorted(self.analytic_Uniprots), sort_keys=True)).hexdigest()
        payload = {'UP_hash' : md5,
                        'sys_hash' : self._MD5hash(),
                        'size' : len(self.analytic_Uniprots),
                        'UPs' : pickle.dumps(self.analytic_Uniprots),
                        'currents' : pickle.dumps((self.current_accumulator, self.node_current)),
                        'voltages' : pickle.dumps(self.UP2circ_and_voltage)}
        dump_object(Dumps.GO_Analysis_memoized, payload)


    def dump_Indep_Linset(self):
        dump_object(Dumps.GO_Indep_Linset, self.Indep_Lapl)


    def undump_Indep_Linset(self):
        self.Indep_Lapl = undump_object(Dumps.GO_Indep_Linset)


    def undump_memoized(self):
        """
        :return: undumped memoized analysis
        :rtype: dict
        """
        return undump_object(Dumps.GO_Analysis_memoized)


    def store(self):
        self.dump_statics()
        self.dump_core()
        self.dump_matrices()
        self.dump_informativities()
        self.dump_inflated_elements()


    def rebuild(self):
        self.get_GO_access()
        self.get_GO_structure()
        self.get_matrixes()
        self.get_GO_Reach()
        if self.ultraspec_cleaned:
            self.filter_out_ultraspecific()
        self.get_Laplacians()
        self.inflate_matrix_and_indexes()


    def load(self):
        """
        Preloads itself from the saved dumps, in case the Filtering system is the same

        """
        Filtr, Initset, corrfactor, ultraspec_cleaned, ultraspec_lvl = self.undump_statics()
        if self.Filtr != Filtr:
            raise Exception("Wrong Filtering attempted to be recovered from storage")
        if self.InitSet != Initset:
            raise Exception("Wrong Initset attempted to be recovered from storage")
        if self.corrfactor != corrfactor:
            raise Exception("Wrong correction factor attempted to be recovered from storage")
        if self.ultraspec_cleaned != ultraspec_cleaned:
            raise Exception("Ultraspecific terms leveling state is not the same in the database as requested")
        if self.ultraspec_lvl != ultraspec_lvl:
            raise Exception("Ultraspecific terms leveling cut-off is not the same in the database as requested")
        self.undump_core()
        self.undump_matrices()
        self.undump_informativities()
        self.undump_inflated_elements()


    def get_GO_access(self):
        """
        Loads all of the relations between the UNIPROTs and GOs as one giant dictionary

        """
        UPs_without_GO = 0
        for UP_DB_ID in self.InitSet:
            UP_Specific_GOs = []
            Root = DatabaseGraph.UNIPORT.get(UP_DB_ID)
            self.UP_Names[UP_DB_ID] = [Root.ID, Root.displayName]
            Node_gen = Root.bothV("is_go_annotation")
            if Node_gen:
                for GO in Node_gen:
                    if GO.Namespace in self.Filtr:
                        GOID = str(GO).split('/')[-1][:-1]
                        UP_Specific_GOs.append(GOID)
                        self.GO2UP[GOID].append(UP_DB_ID)
                        self.SeedSet.add(GOID)
            if not UP_Specific_GOs:
                UPs_without_GO += 1
                print "Warning: UP without GO has been found. Database UP_DB_ID: %s, \t name: %s!!!!!!" % (UP_DB_ID, MG.ID2displayName[UP_DB_ID])
                self.UPs_without_GO.add(UP_DB_ID)
            else:
                self.UP2GO_Dict[UP_DB_ID] = copy.copy(UP_Specific_GOs)

        print 'total number of UPs without a GO annotation: %s out of %s' % (UPs_without_GO, len(self.InitSet))


    def get_GO_structure(self):
        """
        Loads all of the relations between the GOs that are generalisation of the seedList GOs and that are withing the types specified in Filtr

        """
        VisitedSet = set()
        seedList = copy.copy(list(self.SeedSet))
        while seedList:
            ID = seedList.pop()
            VisitedSet.add(ID)
            Local_UpList = []
            Local_Regulation_List = []
            Local_InReg_List = []
            Local_DownList = []
            GONode = DatabaseGraph.GOTerm.get(ID)
            self.GO_Names[ID] = str(GONode.displayName)
            self.GO_Legacy_IDs[ID] = str(GONode.ID)
            self.rev_GO_IDs[str(GONode.ID)] = ID
            for Typ in chain(GOUpTypes, GORegTypes):
                generator = GONode.outV(Typ)
                if not generator:
                    continue  # skip in case GO Node has no outgoing relations to other GO nodes
                for elt in generator:
                    if not elt.Namespace in self.Filtr:
                        continue  # skip in case other GO nodes are of bad type (normally skopes are well-separated, but who knows)
                    subID = str(elt).split('/')[-1][:-1]
                    if subID not in VisitedSet:
                        seedList.append(subID)
                    if Typ in GOUpTypes:
                        Local_UpList.append(subID)
                    else:
                        Local_Regulation_List.append(subID)
                rev_generator = GONode.inV(Typ)
                if not rev_generator:
                    continue
                for elt in rev_generator:
                    if not elt.Namespace in self.Filtr:
                        continue
                    subID = str(elt).split('/')[-1][:-1]
                    if Typ in GOUpTypes:
                        Local_DownList.append(subID)
                    else:
                        Local_InReg_List.append(subID)
            self.Reachable_nodes_dict[ID] = (list(set(Local_UpList)), list(set(Local_Regulation_List)),
                                             list(set(Local_DownList)), list(set(Local_InReg_List)))

        self.All_GOs = list(VisitedSet)
        self.Num2GO = dict( (i, val) for i, val in enumerate(self.All_GOs) )
        self.GO2Num = dict( (val, i) for i, val in enumerate(self.All_GOs) )


    def get_matrixes(self, include_reg = True):
        """
        Builds Undirected and directed adjacency matrices for the GO set and

        :param include_reg: if True, the regulation set is included into the matrix
        :warning: if the parameter above is set to False, get_GO_reach module will be unable to function.
        """

        def build_adjacency():
            """
            Builds undirected adjacency matrix for the GO transitions

            """
            baseMatrix = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
            for node, package in self.Reachable_nodes_dict.iteritems():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.GO2Num[node], self.GO2Num[node2])
                    baseMatrix[idx] = 1
                    idx = (idx[1], idx[0])
                    baseMatrix[idx] = 1

            self.Adjacency_matrix = copy.copy(baseMatrix)



        def build_dir_adj():
            """
            Builds directed adjacency matrix for the GO transitions

            """
            baseMatrix = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
            for node, package in self.Reachable_nodes_dict.iteritems():
                fw_nodes = package[0]
                if include_reg:
                    fw_nodes += package[1]
                for node2 in fw_nodes:
                    idx = (self.GO2Num[node], self.GO2Num[node2])
                    baseMatrix[idx] = 1

            self.dir_adj_matrix = copy.copy(baseMatrix)

        build_adjacency()
        build_dir_adj()


    def _infcalc(self, number):
        """
        returns an entropy given by a number of equiprobable events, where event is the number.

        :param number:
        """
        if number < 1.0:
            raise Exception("Wrong value provided for entropy computation")
        if not self.total_Entropy:
            self.total_Entropy = -log(1/float(len(self.UP2GO_Dict.keys())), 2)
        if number == 1.0:
             return 2*self.total_Entropy
        return pow(-self.corrfactor[0]*self.total_Entropy/log(1/float(number),2), self.corrfactor[1])


    def get_GO_Reach(self):
        """
        Recovers by how many different uniprots each GO term is reached, both in distance-agnostic and distance-specific
        terms.

        """

        def verify_equivalence_of_reaches(step_reach, reach):
            """

            :param step_reach:
            :param reach:
            :raise Exception:
            """
            lendict = {key : [len(val), len(step_reach[key].keys())] for key, val in reach.iteritems()}
            for key, val in lendict.iteritems():
                if val[1] != val[0]:
                    raise Exception('Reach exploration results not equivalent! Please report the error.')


        def special_sum(Dico, filter_funct = lambda x : x+1.0 ):
            """
            Special sum used for the computation of staged informativity of different terms

            :param Dico:
            :param filter_funct:
            :raise Exception:
            """
            summer = 0
            for key, val_list in Dico.iteritems():
                summer += filter_funct(key)*len(val_list)
            return summer


        dir_reg_path = shortest_path(self.dir_adj_matrix, directed = True, method = 'D' )
        dir_reg_path[np.isinf(dir_reg_path)] = 0.0
        dir_reg_path = lil_matrix(dir_reg_path)

        self.GO2UP_Reachable_nodes = dict((el, []) for el in self.Reachable_nodes_dict.keys())
        self.GO2UP_Reachable_nodes.update(self.GO2UP)

        pre_GO2UP_step_reachable_nodes = dict((key, dict((v,0) for v in val)) for key, val in self.GO2UP_Reachable_nodes.iteritems())
        # when called on possibly unencoutenred items, anticipate a default falue of 10 000

        # Now just scan vertical columns and add UP terms attached
        for idx1, idx2 in zip(list(dir_reg_path.nonzero()[0]),list(dir_reg_path.nonzero()[1])):
            self.GO2UP_Reachable_nodes[self.Num2GO[idx2]] += self.GO2UP_Reachable_nodes[self.Num2GO[idx1]]
            if dir_reg_path[idx1, idx2] < 1.0:
                raise Exception("null in non-null patch")
            step_reach_upgrade = dict( (key, val + dir_reg_path[idx1, idx2]) for key, val in pre_GO2UP_step_reachable_nodes[self.Num2GO[idx1]].iteritems())
            for k, v in step_reach_upgrade.iteritems():
                pre_GO2UP_step_reachable_nodes[self.Num2GO[idx2]][k] = min(pre_GO2UP_step_reachable_nodes[self.Num2GO[idx2]].setdefault(k,100000), v)

        for key, val in self.GO2UP_Reachable_nodes.iteritems():
            self.GO2UP_Reachable_nodes[key] = list(set(val))

        verify_equivalence_of_reaches(pre_GO2UP_step_reachable_nodes, self.GO2UP_Reachable_nodes)

        # Now we need to invert the reach to get the set of all the primary and derived GO terms that describe a UP
        self.UP2GO_Reachable_nodes = dict((key, []) for key in self.UP2GO_Dict.keys())
        self.UP2GO_step_Reachable_nodes= dict((key, defaultdict(list)) for key in self.UP2GO_Dict.keys())
        self.GO2UP_step_Reachable_nodes = dict((key, defaultdict(list)) for key in pre_GO2UP_step_reachable_nodes.keys())
        for key, val_dict in pre_GO2UP_step_reachable_nodes.iteritems():
            for k, v in val_dict.iteritems():
                self.GO2UP_step_Reachable_nodes[key][v].append(k)
                self.UP2GO_step_Reachable_nodes[k][v].append(key)
                self.UP2GO_Reachable_nodes[k].append(key)

        # and finally we compute the pure and weighted informativities for each term
        self.GO2_Pure_Inf = dict( (key, self._infcalc(len(val))) for key, val in self.GO2UP_Reachable_nodes.iteritems())
        self.GO2_Weighted_Ent = dict( (key, self._infcalc(special_sum(val_dict))) for key, val_dict in self.GO2UP_step_Reachable_nodes.iteritems())


    def get_Laplacians(self):
        """
        Recovers the Laplacian (information conductance) matrixes for the GO annotation terms.
        For weighted laplacian, currently implements a Max-Ent with custom factor as transition price.

        :warning: for this method to function, get_GO reach function must be run first.
        :warning: accounting for regulatory relation relation between the GO terms is performed if has been done in the adjunction matrix computation

        """
        baseMatrix = -copy.copy(self.dir_adj_matrix)
        nz_list = copy.copy(zip(list(baseMatrix.nonzero()[0]), list(baseMatrix.nonzero()[1])))

        for idx1, idx2 in nz_list:
            minInf = min(self.GO2_Pure_Inf[self.Num2GO[idx1]],self.GO2_Pure_Inf[self.Num2GO[idx2]])
            baseMatrix[idx1, idx2] = -minInf
            baseMatrix[idx2, idx1] = -minInf
            baseMatrix[idx2, idx2] += minInf
            baseMatrix[idx1, idx1] += minInf

        self.Laplacian_matrix = baseMatrix


    def compute_UniprotDict(self):
        """
        Computes the uniprot method requried by some other dictionary

        :return:
        """
        UniprotDict = {}
        for elt in MG.Uniprots:
            node = DatabaseGraph.UNIPORT.get(elt)
            altID = node.ID
            UniprotDict[altID] = (elt, MG.ID2displayName[elt])
            UniprotDict[elt] = altID
        pickle.dump(UniprotDict, file(Dumps.Up_dict_dump,'w'))
        return UniprotDict


    def filter_out_ultraspecific(self):
        """
        Filters out GO terms that are too specific and builds a directed, undirected adjacency maps and laplacian.

        """
        rep_val = self._infcalc(self.ultraspec_lvl)
        self.ultraspec_cleaned = True
        ultraspec_GOs = list( GO for GO, reach in self.GO2UP_Reachable_nodes.iteritems() if len(reach) < self.ultraspec_lvl)
        for GO in ultraspec_GOs:
            self.GO2_Pure_Inf[GO] = rep_val


    def _MD5hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been correctly defined before dump/retrieval
        """
        sorted_initset = sorted(self.InitSet, key=str.lower)
        data = [self.Filtr, sorted_initset, self.corrfactor, self.ultraspec_cleaned, self.ultraspec_lvl]
        md5 = hashlib.md5(json.dumps(data,sort_keys=True)).hexdigest()
        return str(md5)


    def inflate_matrix_and_indexes(self):
        """
        Performs the laplacian matrix inflation to incorporate the uniprots on which we will be running the
        """

        # branching distribution: at least 10x the biggest conductivity of the system, unless ultrasec, in which case ~ ultraspec level
        self.binding_intesity = 10*self._infcalc(self.ultraspec_lvl)
        fixed_Index = self.Laplacian_matrix.shape[0]
        self_connectable_UPs = list(set(self.InitSet) - set(self.UPs_without_GO))
        UP2IDxs = dict((UP, fixed_Index+Idx) for Idx, UP in enumerate(self_connectable_UPs))
        IDx2UPs = dict((Idx, UP) for UP, Idx in UP2IDxs.iteritems())
        self.inflated_Laplacian = lil_matrix((self.Laplacian_matrix.shape[0]+len(self_connectable_UPs),
                                       self.Laplacian_matrix.shape[1]+len(self_connectable_UPs)))
        self.inflated_Laplacian[:self.Laplacian_matrix.shape[0], :self.Laplacian_matrix.shape[1]] = self.Laplacian_matrix

        for UP in self_connectable_UPs:
            for GO in self.UP2GO_Dict[UP]:
                self.inflated_Laplacian[UP2IDxs[UP], UP2IDxs[UP]] += self.binding_intesity
                self.inflated_Laplacian[self.GO2Num[GO], self.GO2Num[GO]] += self.binding_intesity
                self.inflated_Laplacian[self.GO2Num[GO], UP2IDxs[UP]] -= self.binding_intesity
                self.inflated_Laplacian[UP2IDxs[UP], self.GO2Num[GO]] -= self.binding_intesity


        self.inflated_lbl2idx = copy.copy(self.GO2Num)
        self.inflated_lbl2idx.update(UP2IDxs)
        self.inflated_idx2lbl = copy.copy(self.Num2GO)
        self.inflated_idx2lbl.update(IDx2UPs)


    def set_Uniprot_source(self, Uniprots):
        if not set(Uniprots) <= set(self.UP2GO_Dict.keys()):
            ex_pload = 'Following Uniprots either were not in the construction set or have no GOs attached: \n %s' % (set(Uniprots) - set(self.UP2GO_Dict.keys()))
            raise Warning(ex_pload)
        self.analytic_Uniprots = Uniprots


    def build_extended_conduction_system(self, memoized=True, sourced=False, incremental=False, cancellation = True):
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier knowledge flow analysis

        :return: adjusted conduction system
        """
        if not incremental or self.current_accumulator == np.zeros((2,2)):
            self.current_accumulator = lil_matrix(self.inflated_Laplacian.shape)
            self.UP2UP_voltages = {}
            if not sourced:
                self.UP2circ_and_voltage = {}

        for UP1, UP2 in combinations(self.analytic_Uniprots, 2):

            if sourced:
                self.current_accumulator = self.current_accumulator +\
                                           CR.sparse_abs(self.UP2circ_and_voltage[tuple(sorted((UP1, UP2)))][1])
                continue

            Idx1, Idx2 = (self.inflated_lbl2idx[UP1], self.inflated_lbl2idx[UP2])
            pre_reach = self.UP2GO_Reachable_nodes[UP1] + self.UP2GO_Reachable_nodes[UP2] + [UP1] + [UP2]
            reach = [self.inflated_lbl2idx[label] for label in pre_reach]
            current_upper, voltage_diff = CR.get_current_with_reach_limitations(inflated_laplacian = self.inflated_Laplacian,
                                                            Idx_pair = (Idx1, Idx2),
                                                            reach_limiter = reach)
            self.current_accumulator = self.current_accumulator + CR.sparse_abs(current_upper)

            self.UP2UP_voltages[(UP1, UP2)] = voltage_diff

            if memoized:
                self.UP2circ_and_voltage[tuple(sorted((UP1, UP2)))] = (voltage_diff, current_upper)

        if cancellation:
            ln = len(self.analytic_Uniprots)
            self.current_accumulator = self.current_accumulator/(ln*(ln-1)/2)

        if memoized:
            self.dump_memoized()

        index_current = CR.get_current_through_nodes(self.current_accumulator)
        self.node_current = dict( (self.inflated_idx2lbl[idx], val) for idx, val in enumerate(index_current))


    def compute_conduction_system(self, node_current, limit = 0.01):
        charDict = {}
        for GO in self.GO2Num.iterkeys():
            if node_current[GO] > limit:
                charDict[GO] = [ node_current[GO],
                                 self.GO2_Pure_Inf[GO],
                                 len(self.GO2UP_Reachable_nodes[GO])]
        return charDict


    def export_conduction_system(self):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format and flushes it into a file that
        can be viewed with Gephi

        :raise Warning:
        """
        nodecharnames = ['Current', 'Type', 'Legacy_ID', 'Names', 'Pure_informativity', 'Confusion_potential']
        nodechartypes = ['DOUBLE', 'VARCHAR', 'VARCHAR', 'VARCHAR', 'DOUBLE', 'DOUBLE']
        charDict = {}

        for GO in self.GO2Num.iterkeys():
            charDict[GO] = [ str(self.node_current[GO]),
                             'GO', self.GO_Legacy_IDs[GO],
                             self.GO_Names[GO].replace(',','-'),
                             str(self.GO2_Pure_Inf[GO]),
                             str(len(self.GO2UP_Reachable_nodes[GO]))]

        for UP in self.analytic_Uniprots:
            charDict[UP] = [ str(self.node_current[UP]),
                             'UP', self.UP_Names[UP][0],
                             str(self.UP_Names[UP][1]).replace(',','-'),
                             str(self.binding_intesity),
                             '1']

        GDF_exporter = GDF_export_Interface(target_fname = Outputs.GO_GDF_output, field_names = nodecharnames,
                                            field_types = nodechartypes, node_properties_dict = charDict,
                                            mincurrent = 0.01, Idx2Label = self.inflated_idx2lbl, Label2Idx = self.inflated_lbl2idx,
                                            current_Matrix = self.current_accumulator)
        GDF_exporter.write()


    def export_subsystem(self, UP_system, UP_subsystem):
        current_recombinator = self.undump_memoized()
        if not set(UP_system) == set(pickle.loads(current_recombinator['UPs'])):
            raise Exception('Wrong UP system re-analyzed')
        self.UP2circ_and_voltage = pickle.loads(current_recombinator['voltages'])
        self.set_Uniprot_source(UP_subsystem)
        self.build_extended_conduction_system(memoized = False, sourced = True)
        self.export_conduction_system()


    def randomly_sample(self, samples_size, samples_each_size):
        """
        Randomly samples the set

        :param samples_size:
        :param samples_each_size:
        :raise Exception:
        """
        if not len(samples_size) == len(samples_each_size):
            raise Exception('Not the same list sizes!')

        self_connectable_UPs = list(set(self.InitSet) - set(self.UPs_without_GO))

        for sample_size, iterations in zip(samples_size, samples_each_size):
            for i in range(0, iterations):
                shuffle(self_connectable_UPs)
                analytics_UP_list = self_connectable_UPs[:sample_size]
                self.set_Uniprot_source(analytics_UP_list)
                self.build_extended_conduction_system(memoized=False, sourced=False)

                md5 = hashlib.md5(json.dumps( sorted(analytics_UP_list), sort_keys=True)).hexdigest()

                UP_rand_samp.insert({'UP_hash' : md5,
                                     'sys_hash' : self._MD5hash(),
                                     'size' : sample_size,
                                     'UPs' : pickle.dumps(analytics_UP_list),
                                     'currents' : pickle.dumps((self.current_accumulator, self.node_current)),
                                     'voltages' : pickle.dumps(self.UP2UP_voltages)})

                print 'Random ID: %s \t Sample size: %s \t iteration: %s\t compops: %s \t time: %s ' %(self.r_ID,
                        sample_size, i, "{0:.2f}".format(sample_size**2/2/self._time()), self.pretty_time())


    def get_indep_linear_groups(self):
        self.Indep_Lapl = lil_matrix((len(self.All_GOs), len(self.All_GOs)))
        for GO_list in self.UP2GO_Reachable_nodes.itervalues():
            for GO1, GO2 in combinations(GO_list,2):
                idx1, idx2 = (self.GO2Num[GO1], self.GO2Num[GO2])
                self.Indep_Lapl [idx1, idx2] += -1
                self.Indep_Lapl [idx2, idx1] += -1
                self.Indep_Lapl [idx2, idx2] += 1
                self.Indep_Lapl [idx1, idx1] += 1


if __name__ == '__main__':
    filtr = ['biological_process']

    KG = GO_Interface(filtr, MG.Uniprots, (0.3, 2), True, 3)
    KG.rebuild()
    print KG.pretty_time()
    KG.store()
    print KG.pretty_time()
    # experimental = ['881579','65094', '925081', '500332', '915530', '456374']
    #
    # KG.load()
    # print KG.pretty_time()
    # KG.get_indep_linear_groups()
    # KG.dump_Indep_Linset()
    #
    #
    #
    #
    # KG.randomly_sample([10], [5])
    #
    # KG.set_Uniprot_source(experimental)
    # KG.build_extended_conduction_system()
    # KG.export_conduction_system()
    #
    # KG.export_subsystem(experimental, ['881579','65094'])




    # TODO: build a theoretic "surprise" level estimator, comapare with experimental "surprise" of the GO terms with certain confusion/informativity observation

    # TODO: estimate how likely it is to get a term with informativity x at a size N and what is the expectend number of terms, use the Benjamin-Hoechberg method
            # to remove the most unlikely ones.

    # Non-trivial interesting GOs: reach between 3 and 200. For them, we should calculate the hidden strong importance
    # terms, i.e. routing over X percent of information => UP importance for GO terms.

    # loading takes 1-6 seconds.
    # fill for reach only is done in 2 seconds,
    # tepping takes another 15,
    # inverting + info computation - 1 more second
    # Laplacian building =>

    # full computation - 3 minutes 18 seconds; save 7 seconds, retrieval - 3 seconds

    # data_array = np.array([log(val) for val in KG.GO2_Pure_Inf.itervalues()])
    # hist(data_array, 100, log=True)
    # show()

    # UP confusion matrix by GO Term?
        # => coefficient of corelation?
    # GO specificity (i.e. how fully the given UP set covers a given GO term )
    # GO recall (i.e. how many UPs a given GO term covers from a given set? )
    # GO confusion power: i.e. knowing just that GO terms, how many terms out of this set would we retrieve? in case of
    # an equiprobable pulling out of all the UP terms you know are annotated with such a GO?

    # How can we scale it with introducing credibility values for the Uniprots?

    # Sorting on the probability of confusion with a random sample, of addition of a random sample of a given size to the
    # GO analysis system?

    # what would be the probability of a set with a given informativity popping out on random for a given selection of terms?

    # in order to perform the correct implementation, we need to draw a graph of the correlation between the confusion
    # of the minimal informativity GO term of a set and the tension required to transmit the information

    # Problem: ultraspecificity requires a pretty strong re-manipulation of both GOs and matrix construction routines after
    # the matrices have been build.

    # Wouldn't it be easier just to repress the informativity they induce by modifying the correction factors in the
    # information functions?
        # => Now, the distortion the uinduce is more then 10 times stronger. In addition the presence of this terms
        # would greatly slow the matrix computation performance. (there are almost 3000 of them: a third of all the GO terms)
    # What we could do is just to ceil the informativity of terms on the x-term informaitivity
