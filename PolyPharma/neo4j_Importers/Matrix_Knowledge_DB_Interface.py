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
from PolyPharma.configs import Dumps
from PolyPharma.neo4j_analyzer.knowledge_access import acceleratedInsert
from PolyPharma.neo4j_analyzer.IO_Routines import dump_object, write_to_csv, undump_object

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
                                                                len(matrix.nonzero()[0]), np.min(matrix), np.max(matrix))


class GO_Interface(object):
    """
    General calss to recover all the informations associated with GO from database and buffer them for further use.
    """

    def __init__(self, Filter, Uniprot_Node_IDs, corrfactor):
        self.Filtr = Filter
        self.InitSet = Uniprot_Node_IDs
        self.corrfactor = corrfactor
        self.init_time = time()
        self.partial_time = time()

        self.UP2GO_Dict = {}
        self.GO2UP = defaultdict(list)
        self.SeedSet = set()
        self.All_GOs = []
        self.GO2Num = {}
        self.Num2GO = {}
        self.total_Entropy = None

        self.UP2GO_Reachable_nodes = {}
        self.GO2UP_Reachable_nodes = {}
        self.UP2GO_step_Reachable_nodes = {}
        self.GO2UP_step_Reachable_nodes = {}
        self.GO2_Pure_Inf = {}
        self.GO2_Weighted_Ent = {}

        self.GO_Names = {}
        self.GO_Legacy_IDs = {}
        self.rev_GO_IDs = {}

        self.Adjacency_matrix = np.zeros((2, 2))
        self.dir_adj_matrix = np.zeros((2, 2))
        self.Laplacian_matrix = np.zeros((2, 2))
        self.Weighted_Laplacian_matrix = np.zeros((2, 2))
        self.Sign_retaining_matrix = np.zeros((2, 2))

        self.TimesReached = {}
        self.accelerationDict = {}
        self.Reverse_Dict = {}
        self.GO_Node_ID2Reach = {}


    def time(self):
        """
        Times the execution

        :return: tuple containing the time since the creation of the Matrix_getter object and since the last cal of function formatted as string
        :rtype: str
        """
        it, pt = (round(time() - self.init_time), round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (int(it) / 60, it % 60, int(pt) / 60, pt % 60)
        self.partial_time = time()
        return pload


    def dump_statics(self):
        dump_object(Dumps.GO_builder_stat, (self.Filtr, self.InitSet ,self.corrfactor))


    def undump_statics(self):
        return undump_object(Dumps.GO_builder_stat)


    def dump_core(self):
        dump_object(Dumps.GO_dump, (self.UP2GO_Dict, self.GO2UP, self.SeedSet, self.Reachable_nodes_dict, self.GO_Names,
                                    self.GO_Legacy_IDs, self.rev_GO_IDs, self.All_GOs, self.GO2Num, self.Num2GO))


    def undump_core(self):
        self.UP2GO_Dict, self.GO2UP, self.SeedSet, self.Reachable_nodes_dict, self.GO_Names, self.GO_Legacy_IDs,\
        self.rev_GO_IDs, self.All_GOs, self.GO2Num, self.Num2GO = undump_object(Dumps.GO_dump)


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


    def store(self):
        self.dump_statics()
        self.dump_core()
        self.dump_matrices()
        self.dump_informativities()


    def rebuild(self):
        self.get_GO_access()
        self.get_GO_structure()
        self.get_matrixes()
        self.get_GO_Reach()


    def load(self):
        """
        Preloads itself from the saved dumps, in case the Filtering system is the same

        """
        Filtr, Initset, corrfactor = self.undump_statics()
        if self.Filtr != Filtr:
            raise Exception("Wrong Filtering attempted to be recovered from storage")
        if self.InitSet != Initset:
            raise Exception("Wrong Initset attempted to be recovered from storage")
        if self.corrfactor != corrfactor:
            raise Exception("Wrong correction factor attempted to be recovered from storage")
        self.undump_core()
        self.undump_matrices()

    def get_GO_access(self):
        """
        Loads all of the relations between the UNIPROTs and GOs as one giant dictionary

        """
        UPs_without_GO = 0
        for UP_DB_ID in self.InitSet:
            UP_Specific_GOs = []
            Root = DatabaseGraph.UNIPORT.get(UP_DB_ID)
            Node_gen = Root.bothV("is_go_annotation")
            if Node_gen:
                for GO in Node_gen:
                    if GO.Namespace in self.Filtr:
                        GOID = str(GO).split('/')[-1][:-1]
                        UP_Specific_GOs.append(GOID)
                        self.GO2UP[GOID].append(UP_DB_ID)
                        self.SeedSet.add(GOID)
            if UP_Specific_GOs == []:
                UPs_without_GO += 1
                print "Warning: UP without GO has been found. Database UP_DB_ID: %s, \t name: %s!!!!!!" % (UP_DB_ID, MG.ID2displayName[UP_DB_ID])
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


    def equient(self, number):
        """
        returns an entropy given by a number of equiprobable events, where event is the number.

        :param number:
        """
        if number < 1.0:
            raise Exception("Wrong value provided for entropy computation")
        if not self.total_Entropy:
            self.total_Entropy = -log(1/float(len(self.UP2GO_Dict.keys())),2)
        return pow(-self.total_Entropy/log(1/float(number),2), self.corrfactor)


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
        print dir_reg_path.nonzero()

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
        self.GO2_Pure_Inf = dict( (key, self.equient(len(val))) for key, val in self.GO2UP_Reachable_nodes.iteritems())
        self.GO2_Weighted_Ent = dict( (key, self.equient(special_sum(val_dict))) for key, val_dict in self.GO2UP_step_Reachable_nodes.iteritems())


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
            baseMatrix[idx1, idx2] = -minInf # TODO: figrue out if entropic information computation should be done here or in a special function.
            baseMatrix[idx2, idx1] = -minInf # TODO: reduction of information to plafonate it at the triad-level
            baseMatrix[idx2, idx2] += minInf
            baseMatrix[idx1, idx1] += minInf

        self.Laplacian_matrix = baseMatrix


    def get_GO_Informativities(self):
        """
        Here calculated without any information on regulation
        ..todo: introduce the iformation computation with regulation,. Outiline: transmit to each higher-level node the
                number of nodes attainable from each node. I.e. instead of an exploration rooted at each node, use
                a simultaneous exploration of a whole tree.

        """
        i = 0
        l = len(self.UP2GO_Dict)
        for key in self.UP2GO_Dict.keys():
            i += 1
            print 'entering', float(i)/float(l), self.time()
            toVisit = copy.copy(self.UP2GO_Dict[key])
            visited = []
            while toVisit:
                elt = toVisit.pop()
                vs = acceleratedInsert(self.UP2GO_Dict, self.accelerationDict, elt)
                visited += vs
            visited = list(set(visited))
            for elt in visited:
                if elt not in self.TimesReached.keys():
                    self.Reverse_Dict[elt] = []
                    self.TimesReached[elt] = 0
                self.Reverse_Dict[elt].append(key)
                self.TimesReached[elt] += 1
        # Fle=file('GO_Informativities.dump','w')
        # pickle.dump(TimesReached,Fle)  #TODO": correct dumping here
        # Fle2=file('accDict.dump','w')
        # pickle.dump(accelerationDict,Fle2)  #TODO": correct dumping here
        # Fle3=file('Reverse_dict.dump','w')
        # pickle.dump(Reverse_Dict, Fle3)  #TODO": correct dumping here


    def compute_UniprotDict(self):
        """


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


if __name__ == '__main__':
    filtr = ['biological_process']

    KG = GO_Interface(filtr, MG.Uniprots, 1)
    KG.rebuild()
    KG.store()
    print KG.time()

    KG.load()
    print KG.time()
    KG.get_GO_Reach()
    print KG.time()

    # Non-trivial interesting GOs: reach between 3 and 200. For them, we should calculate the hidden strong importance
    # terms, i.e. routing over X percent of information => UP importance for GO terms.

    # loading takes 1 second.
    # fill for reach only is done in 2 seconds,
    # tepping takes another 15,
    # inverting + info computation - 1 more second

    data_array = np.array(list(KG.GO2_Pure_Inf.values()))
    hist(data_array, 100, log=True)
    show()