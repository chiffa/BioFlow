__author__ = 'ank'
"""
This module contains all the routines that are respojnsible for pulling
the matrixes out of the neo4j graph and processing them.
"""

import os
import hashlib
import json
import itertools
import pickle
import string
import numpy as np
from copy import copy
from time import time
from random import shuffle, sample
from collections import defaultdict
from pprint import PrettyPrinter
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components

from BioFlow.configs2 import Chromosome_source, Chromosome_file_filter
from BioFlow.configs2 import edge_type_filters, Adjacency_Martix_Dict, Conductance_Matrix_Dict, Dumps, Outputs, Interactome_rand_samp
from BioFlow.Utils.GDF_export import GDF_export_Interface
from BioFlow.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from BioFlow.neo4j_analyzer.DB_IO_Routines import reaction_participant_getter, expand_from_seed, Erase_custom_fields
from BioFlow.neo4j_analyzer.IO_Routines import write_to_csv, dump_object, undump_object
from BioFlow.neo4j_analyzer import Conduction_routines as CR

# Debug Log:
#       Main connex set is 24k nodes, with 4331 UP links and 1051 Hint links
#       With the full UP import, the main connex set is 25k Nodes, with 4336 UP Links and 1191 HiNT links
#       With forbidding overloaded items: 24k771 Nodes, 4293 UP Links, 1186 HiNT links


class MatrixGetter(object):
    """
    Builds an object that is used ans an interface between matrix representation of an interactome and the

    :param Connexity_Aware: if set to True, loads only the elements of the main connex set of the database
    :type Connexity_Aware: bool
    :param full_impact: if set to True, will compare the names of entities and link the entities that has the same \
            name with the "Possibily the same" relation
    :type full_impact: bool
    """

    # List of all the reaction types present in the DatabaseGraph that will be
    # used as roots to build the interaction network (not all nodes are necessary
    # within the connex part of the graph)
    ReactionsList = [DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]


    def __init__(self, Connexity_Aware, full_impact):
        self.Connexity_Aware = Connexity_Aware
        self.full_impact = full_impact

        self.init_time = time()
        self.partial_time = time()

        self.Ajacency_Matrix = np.zeros((4,4))
        self.Conductance_Matrix = np.zeros((4,4))    # This is just non-normalized laplacian matrix

        self.adj_eigenvects = np.zeros((4,4))
        self.adj_eigenvals = np.zeros((4,4))
        self.cond_eigenvects = np.zeros((4,4))
        self.cond_eigenvals = np.zeros((4,4))

        self.NodeID2MatrixNumber = {}
        self.MatrixNumber2NodeID = {}
        self.ID2displayName = {}
        self.ID2LegacyId = {}
        self.ID2Type = {}
        self.ID2Localization = {}
        self.Uniprots = []
        self.Uniprot_complete = []
        self.Uniprot_attachments = {} # currently maintained for legacy reasons
        self.Uniprot_Mat_idxs = []

        self.ReactLinks = []
        self.InitSet = []
        self.GroupLinks = {}
        self.GroupSet = []
        self.SecLinks = {}
        self.SecSet = []
        self.UP_Links = {}
        self.UPSet = []
        self.HiNT_Links = {}
        self.FullSet = []
        self.Super_Links = {}
        self.ExpSet = []

        self.Matrix_reality = Dumps.matrix_corrs

        char_set = string.ascii_uppercase + string.digits
        self.r_ID = ''.join(sample(char_set*6, 6))

        self.analytic_Uniprots = []
        self.UP2UP_voltages = {}
        self.UP2voltage_and_circ = {}
        self.current_accumulator = np.zeros((2,2))
        self.node_current = {}

        self.UP2Chrom = {}
        self.Chrom2UP = defaultdict(list)

        self.incomplete_compute = False  # used in case of sparse sampling
        self.background = None

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


    def dump_Matrices(self):
        """
        dumps self.Ajacency_Matrix and self.Conductance_Matrix
        """
        dump_object(Dumps.ValMat, self.Ajacency_Matrix)
        dump_object(Dumps.ConMat, self.Conductance_Matrix)


    def undump_Matrices(self):
        """
        undumps self.Ajacency_Matrix and self.Conductance_Matrix
        """
        self.Ajacency_Matrix = undump_object(Dumps.ValMat)
        self.Conductance_Matrix = undump_object(Dumps.ConMat)


    def dump_Eigens(self):
        """
        dumps self.adj_eigenvals and self.Conductance_Matrix and writes them to csv
        """
        write_to_csv(Dumps.eigen_VaMat, self.adj_eigenvals)
        write_to_csv(Dumps.eigen_ConMat, self.cond_eigenvals)
        dump_object(Dumps.val_eigen, (self.adj_eigenvals, self.adj_eigenvects))
        dump_object(Dumps.cond_eigen, (self.cond_eigenvals, self.cond_eigenvects))


    def undump_Eigens(self):
        """
        undumps self.adj_eigenvals and self.Conductance_Matrix
        """
        self.adj_eigenvals, self.adj_eigenvects = undump_object(Dumps.val_eigen)
        self.cond_eigenvals, self.cond_eigenvects = undump_object(Dumps.cond_eigen)



    def dump_Maps(self):
        """
        dumps all the elements required for the mapping between the types and ids of database entries and matrix columns
        """
        dump_object(Dumps.matrix_corrs,
                         (self.NodeID2MatrixNumber, self.MatrixNumber2NodeID,
                          self.ID2displayName, self.ID2LegacyId, self.ID2Type, self.ID2Localization,
                          self.Uniprots, self.Uniprot_complete, self.Uniprot_attachments, self.UP2Chrom, self.Chrom2UP,
                          self.Uniprot_Mat_idxs))


    def undump_Maps(self):
        """
        undumps all the elements required for the mapping between the types and ids of database entries and matrix columns
        """
        self.NodeID2MatrixNumber, self.MatrixNumber2NodeID, self.ID2displayName, self.ID2LegacyId, self.ID2Type,\
        self.ID2Localization, self.Uniprots, self.Uniprot_complete, self.Uniprot_attachments, self.UP2Chrom,\
        self.Chrom2UP, self.Uniprot_Mat_idxs = undump_object(Dumps.matrix_corrs)



    def dump_memoized(self):
        md5 = hashlib.md5(json.dumps( sorted(self.analytic_Uniprots), sort_keys=True)).hexdigest()
        payload = {'UP_hash' : md5,
                        'sys_hash' : self._MD5hash(),
                        'size' : len(self.analytic_Uniprots),
                        'UPs' : pickle.dumps(self.analytic_Uniprots),
                        'currents' : pickle.dumps((self.current_accumulator, self.node_current)),
                        'voltages' : pickle.dumps(self.UP2voltage_and_circ)}
        dump_object(Dumps.Interactome_Analysis_memoized, payload)


    def undump_memoized(self):
        """
        :return: undumped memoized analysis
        :rtype: dict
        """
        return undump_object(Dumps.Interactome_Analysis_memoized)


    def full_load_LS(self):

        """
        Performs the loading of all the objects and the relations in the database according to the following algorithm:
            - get all the reaction nodes, i.e. reactions specified in the self.Reaction_List.
            - For each reaction, add the participants to the "seeds" list and create a set of groups linking them all together
            - for each element of hte "seeds" list, crawl the objects reachable from this seed by specific relations,
            add those new elements to the list and remember as reachable from the current element of the seed list.

            The relations are crawled in the following order:
                * "Group"
                * "Same"
                * "Contact_interaction" (repeated 5 times to increase as much as possible the main connex size, even for the knowledge groups that are badly represented in the interactome database)
                * "HiNT Contact_interaction"
                * "possibly_same"
            The exact relations in each group of relation abbreviations are specified in the configs.py file,

        This function fills in the self.Set and self.Links elements of the MatrixGetter object
        """

        def get_Reaction_blocks():
            """
            Recovers the blocks if interaction that are due to a common set of reactions
            for the elements. They will be used as roots to build the complete interaction
            tree later on.

            :return:* Reagent Clusters (list of lists, where each list represents a set of
                    reagents participating to the same reaction),
                    * Seeds list,
                    * number of elements in the seeds list
            :rtype: 3- tuple
            """
            ReagentClusters = []
            Seeds = set()
            count = 0
            for ReactionType in self.ReactionsList:
                if not ReactionType.get_all():
                    continue
                for Reaction in ReactionType.get_all():
                    if Reaction is None:
                        continue
                    LocalList, cnt = reaction_participant_getter(Reaction, self.Connexity_Aware)
                    count += cnt

                    if len(LocalList) > 1:
                        ReagentClusters.append(copy(LocalList))
                        Seeds.update(LocalList)

            return ReagentClusters, Seeds, count

        def get_expansion( SubSeed, edge_type_filter):
            '''
            Recovers all the nodes reached from the SubSeed according to a the relations listed
            in edge_type_filter

            :param SubSeed: set of NodeIDs that serve as a root for searching further relations
            :type SubSeed: set of NodeIDs

            :param edge_type_filter: type of relations according to which the graph will be explored
            :type edge_type_filter: relations of the type bulbsGraph.Relationtype (these are well python objects)

            :return Clusters: * Clusters (Dictionary of lists, where key is the NodeID from the SubSeed)
                              * SuperSeed (List of NodeIDs attained from the SubSet by the relation types from the edge_type_filter)
                              * count (number of elements added to the SuperSeed compared to the SubSeed)
            :rtype: 3- tuple
            '''
            Clusters = {}
            SuperSeed = set()
            SuperSeed.update(SubSeed)
            count = 0
            for element in SubSeed:
                LocalList, count_increase = expand_from_seed(element, edge_type_filter, self.Connexity_Aware)
                if len(LocalList) > 0:
                    Clusters[element] = copy(LocalList)
                    SuperSeed.update(LocalList)
                    count += count_increase
            return Clusters, SuperSeed, count


        def characterise(name, Links, Group, count):
            """
            Prints a charectristics of a return object from the get_Reaction_blocks or get_expansion

            :param name: name of the object characterised
            :type name: str
            :param Links: element 1 of the return tuple (Links between objects)
            :param Group: element 2 of the return tuple (SuperSeed)
            :type Group: set or list
            :param count: element 3 of the return tuple (count)
            :type count: int
            """
            print '===========>', name, '<==========='
            print len(Links), len(Group), count
            print self.pretty_time()

        ################################################################################################################

        # TODO: this step is inherently recursive and must be formalized as such, intaking only the database and
        #       data types defined in the configs files

        self.ReactLinks, self.InitSet, c = get_Reaction_blocks()
        characterise('Reactions', self.ReactLinks, self.InitSet, c)

        self.GroupLinks, self.GroupSet, c = get_expansion(self.InitSet, edge_type_filters["Group"])
        characterise('Groups', self.GroupLinks, self.GroupSet, c)

        self.SecLinks, self.SecSet, c = get_expansion(self.GroupSet, edge_type_filters["Contact_interaction"])
        characterise('Secondary Links', self.SecLinks, self.SecSet, c)

        for i in range(0, 5):
            SecLinks2, SecSet2, c = get_expansion(self.SecSet, edge_type_filters["Contact_interaction"])
            self.SecSet = SecSet2
            self.SecLinks = SecLinks2
            characterise('Secondary Links ' + str(i) + ' ', self.SecLinks, self.SecSet, c)

        self.UP_Links, self.UPSet, c = get_expansion(self.SecSet, edge_type_filters["Same"])
        characterise('Uniprot Links', self.UP_Links, self.UPSet, c)

        self.HiNT_Links, self.pre_FullSet, c = get_expansion(self.UPSet, edge_type_filters["HiNT_Contact_interaction"])
        characterise('HiNT Links', self.HiNT_Links, self.pre_FullSet, c)

        for i in range(0, 5):
            HiNT_Links2, FullSet2, c = get_expansion(self.pre_FullSet, edge_type_filters["HiNT_Contact_interaction"])
            self.pre_FullSet = FullSet2
            self.HiNT_Links = HiNT_Links2
            characterise('HiNT Links ' + str(i) + ' ', self.HiNT_Links, self.pre_FullSet, c)

        self.BioGRID_Links, self.FullSet, c = get_expansion(self.pre_FullSet, edge_type_filters["BioGRID_Contact_interaction"])
        characterise('BioGRID Links', self.BioGRID_Links, self.FullSet, c)

        for i in range(0, 5):
            BioGRID_Links2, FullSet2, c = get_expansion(self.FullSet, edge_type_filters["BioGRID_Contact_interaction"])
            self.FullSet = FullSet2
            self.BioGRID_Links = BioGRID_Links2
            characterise('HiNT Links ' + str(i) + ' ', self.BioGRID_Links, self.FullSet, c)

        self.Super_Links, self.ExpSet, c = get_expansion(self.FullSet, edge_type_filters["possibly_same"])
        characterise('Looks_similar Links', self.Super_Links, self.ExpSet, c)


    def map_rows_to_names(self,):
        """ Maps Node Database IDs, Legacy IDs, display names and types to matrix row/column indexes; """

        def request_location(Location_Buffer_Dict, location):
            """
            Just a Buffered lookup of location, since the number of cellular location is relatively small
            (~80), it makes sense to buffer the IOs on it.
            Normally should be moved out as a buffering decorator

            :type Location_Buffer_Dict: dict
            :param Location_Buffer_Dict: Buffered location
            :param location: location Node Lagacy ID we are willing to verify
            :return: displayName of the requested location
            """
            location = str(location)
            if location in Location_Buffer_Dict.keys():
                return Location_Buffer_Dict[location]
            else:
                generator = DatabaseGraph.Location.index.lookup(ID = location)
                if generator is not None:
                    for elt in generator:
                        Location_Buffer_Dict[location] = str(elt.displayName)
                        return str(elt.displayName)

        ######################################################################################################################


        counter = 0
        LocationBufferDict = {}

        self.NodeID2MatrixNumber = {}
        self.MatrixNumber2NodeID = {}

        for ID in self.Highest_Set:
            self.NodeID2MatrixNumber[ID] = counter
            self.MatrixNumber2NodeID[counter] = ID
            Vertex = DatabaseGraph.vertices.get(ID)
            self.ID2displayName[ID] = Vertex.displayName
            self.ID2Type[ID] = Vertex.element_type
            self.ID2LegacyId[ID] = Vertex.ID
            if Vertex.element_type == "UNIPROT":
                self.Uniprots.append(ID)
                self.Uniprot_Mat_idxs.append(counter)
            if Vertex.localization is not None:
                self.ID2Localization[ID] = request_location(LocationBufferDict, Vertex.localization)
            counter += 1

        self.Uniprot_complete += self.Uniprots
        UP_generator = DatabaseGraph.UNIPORT.get_all()
        if UP_generator:
            for UP_Node in UP_generator:
                ID = str(UP_Node).split('/')[-1][:-1]
                if ID not in self.Uniprots:
                    self.Uniprot_complete.append(ID)
                    self.ID2displayName[ID] = UP_Node.displayName
                    self.ID2Type[ID] = UP_Node.element_type
                    self.ID2LegacyId[ID] = UP_Node.ID

        self.Uniprot_complete = list(set(self.Uniprot_complete))



    def fast_row_insert(self, element, index_type):
        """
        Performs an correct insertion of an edge to the matrix.

        :param element: tuple of indexes designating elements we are willing to link
        :type element: 2- tuple of ints
        :param index_type: type of the insert, so that the matrix coefficient can be looked up in the Adjacency_Martix_Dict
                            or Conductance_Martix_Dict from the configs file
        :type index_type: str
        """
        self.Ajacency_Matrix[element[0], element[1]] = min(self.Ajacency_Matrix[element[0], element[1]] + Adjacency_Martix_Dict[index_type], 1)
        self.Ajacency_Matrix[element[1], element[0]] = min(self.Ajacency_Matrix[element[1], element[0]] + Adjacency_Martix_Dict[index_type], 1)

        self.Conductance_Matrix[element[0], element[1]] -= Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[1], element[0]] -= Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[1], element[1]] += Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[0], element[0]] += Conductance_Matrix_Dict[index_type]


    def create_val_matrix(self):
        """
        Creates anew self.Ajacency_Matrix and self.Conductance_Matrix
        """

        self.full_load_LS()

        self.Highest_Set = self.FullSet

        if self.full_impact:
            self.Highest_Set = self.ExpSet

        print "building correspondances", self.pretty_time()

        self.map_rows_to_names()

        print "building the ValMatrix", self.pretty_time()

        loadLen = len(self.Highest_Set)
        self.Ajacency_Matrix = lil_matrix((loadLen, loadLen))
        self.Conductance_Matrix = lil_matrix((loadLen, loadLen))

        for group in self.ReactLinks:
            for elt in itertools.combinations(group, 2):
                element = (self.NodeID2MatrixNumber[elt[0]], self.NodeID2MatrixNumber[elt[1]])
                self.fast_row_insert(element, "Reaction")

        for key in self.GroupLinks.keys():
            for val in self.GroupLinks[key]:
                element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                self.fast_row_insert(element, "Group")

        for key in self.SecLinks.keys():
            for val in self.SecLinks[key]:
                element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                self.fast_row_insert(element, "Contact_interaction")

        for key in self.UP_Links.keys():
            for val in self.UP_Links[key]:
                element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                self.fast_row_insert(element, "Same")

        for key in self.HiNT_Links.keys():
            for val in self.HiNT_Links[key]:
                element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                self.fast_row_insert(element, "Contact_interaction")

        for key in self.BioGRID_Links.keys():
            for val in self.BioGRID_Links[key]:
                element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                self.fast_row_insert(element, "weak_contact")

        if self.full_impact:
            for key in self.Super_Links.keys():
                for val in self.Super_Links[key]:
                    element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                    self.fast_row_insert(element, "possibly_same")


    def get_eigenspectrum(self, numberEigvals):
        """
        Recovers the eigenspectrum associated to the *n* biggest eigenvalues, where *n* is specified by numberEigvals.
            If the Adjacency and conductance matrix haven't been preloaded first, will raise an Exception

        :param numberEigvals: specifies how many biggest eigenvalues we are willing to get.
        :type numberEigvals: int
        :raise Exception: "Matrix must be pre-loaded first" if self.Ajacency_Matrix and self.Conductance_Matrix have not
                            been computed anew or pre-loaded first
        """
        if self.Conductance_Matrix.shape == (4,4):
            raise Exception("Matrix must be pre-loaded first")

        print "entering eigenvect computation", self.pretty_time()

        self.adj_eigenvals, self.adj_eigenvects = eigsh(self.Ajacency_Matrix, numberEigvals)
        self.cond_eigenvals, self.cond_eigenvects = eigsh(self.Conductance_Matrix, numberEigvals)

        print self.adj_eigenvals
        print '<======================>'
        print self.cond_eigenvals
        print '<======================>'
        print np.all(eigsh(self.Conductance_Matrix)[0] > 0)

        print self.pretty_time()


    def get_normalized_laplacian(self, fudge = 1E-18):
        """
        Returns the normalized conductance matrix (Conductance matrix is actually the matrix laplacian of the associated
        Graph)

        :param fudge: Cancells out the action of eigenvalues that are too low.
        :return: normalized matrix
        """
        d, V = np.linalg.eigh(self.Conductance_Matrix)
        D = np.diag(1./np.sqrt(d + fudge))
        return np.dot(np.dot(V, D), V.T)


    def Write_Connexity_Infos(self):
        """
            Writes the infos about connexity of different components of a graph into the database for the future use.
            This execution is the main reason for the existance of the Adjacency Matrix.

            :warning: This process has to ber re-run each time the underlying database is changed in a way that migh affect the main connex graph
        """
        CCOmps = connected_components(self.Ajacency_Matrix, directed = False)

        counters = np.zeros((CCOmps[0], 1))
        for elt in range(0, len(CCOmps[1])):
            counters[CCOmps[1][elt], 0] += 1
        major_Index = np.argmax(counters)

        ln = len(CCOmps[1])
        for i in range(0, len(CCOmps[1])):
            print 'Marking graph main connex elements: %s done' % str("{0:.2f}".format(float(i)/float(ln)*100))
            if CCOmps[1][i] == major_Index:
                Node = DatabaseGraph.vertices.get(self.MatrixNumber2NodeID[i])
                Node.custom = 'Main_Connex'
                Node.main_connex = True
                Node.save()
        print "Marking of %s nodes for connexity was done in %s" % (str(ln), str(self.pretty_time()))


    def full_rebuild(self):
        """
        Performs the initial loading routines that set up in place the sytem of dump files to allow fast loading in the
        future
        """

        CA  = self.Connexity_Aware
        self.Connexity_Aware = False

        self.create_val_matrix()

        Erase_custom_fields()
        self.Write_Connexity_Infos()

        self.compute_Uniprot_Attachments()
        self.Connexity_Aware = CA
        self.create_val_matrix()
        self.get_eigenspectrum(100)

        self.map_UPs_to_chromosomes()

        self.dump_Maps()
        self.dump_Matrices()
        self.dump_Eigens()


    def fast_load(self):
        """
        Preloads mappings, matrices and eigenvalues
        """
        self.undump_Maps()
        self.undump_Matrices()
        self.undump_Eigens()


    def get_descriptor_for_index(self, index):
        """
        Recovers a descriptor set for a given index in the current matrix mapping

        :param index: idenx for which return a descirptor
        :return: Type, displayName and if a localization is given, returns display name too.
        :rtype: tuple
        """
        if self.MatrixNumber2NodeID[index] in self.ID2Localization.keys():
            return (self.ID2Type[self.MatrixNumber2NodeID[index]],
                    self.ID2displayName[self.MatrixNumber2NodeID[index]],
                    self.ID2Localization[self.MatrixNumber2NodeID[index]])
        else:
            return (self.ID2Type[self.MatrixNumber2NodeID[index]],
                    self.ID2displayName[self.MatrixNumber2NodeID[index]])


    def compute_Uniprot_Attachments(self):
        """
        Computes the dictionarry of attachements between the Uniprots and Reactome proteins
        """
        for SP_Node_ID in self.Uniprots:
            Vertex = DatabaseGraph.UNIPORT.get(SP_Node_ID)
            Generator = Vertex.bothV("is_same")
            if Generator is not None:
                self.Uniprot_attachments[SP_Node_ID] = []
                for item in Generator:
                    ID = str(item).split('/')[-1][:-1]
                    self.Uniprot_attachments[SP_Node_ID].append(ID)
                print 'attached %s Reactome proteins to the node %s' %(len(self.Uniprot_attachments[SP_Node_ID]), SP_Node_ID)
            else:
                print 'No attachement for the node %s' % SP_Node_ID


    def hacky_corr(self):
        """
        Hacky method that should remain unused but by the devs.
        Generate the Uniprot_Mat_idxs from the loaded Uniprot List. Prefents a 90-minute full database reload to compute
        Sthing quickly.
        """
        self.undump_Maps()
        self.Uniprot_Mat_idxs = []
        for SP_Id in self.Uniprots:
            Node = DatabaseGraph.UNIPORT.get(SP_Id)
            if Node.main_connex:
                self.Uniprot_Mat_idxs.append(self.NodeID2MatrixNumber[SP_Id])
        print len(self.Uniprot_Mat_idxs)
        self.dump_Maps()


    def _MD5hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been correctly defined before dump/retrieval
        """
        sorted_initset = sorted(self.NodeID2MatrixNumber.keys(), key=str.lower)
        data = [self.Connexity_Aware, sorted_initset, self.full_impact, self.background]
        md5 = hashlib.md5(json.dumps(data, sort_keys=True)).hexdigest()
        return str(md5)


    def set_Uniprot_source(self, Uniprots):
        """
        Sets the Uniprots on which the circulation computation routines will be performed by the otehr methods.
        Avoids passing as argument large lists of parameters.

        :param Uniprots: List of node IDs of the uniprots on which we would like to perform current computations
        :raise Warning: if the uniprots were not present in the set of GOs for which we built the system or had no GO attached to them
        """
        if not set(Uniprots) <= set(self.NodeID2MatrixNumber.keys()):
            ex_pload = 'Following Uniprots were not retrieved upon the circulation matrix construction: \n %s' % (set(Uniprots) - set(self.NodeID2MatrixNumber.keys()))
            print Warning(ex_pload)
        self.analytic_Uniprots =self.analytic_Uniprots = [ uniprot for uniprot in Uniprots if uniprot in self.NodeID2MatrixNumber.keys()]


    def build_extended_conduction_system(self, memoized=True, sourced=False, incremental=False, cancellation=True, sparse_samples=False):
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier knowledge flow analysis

        :param memoized: if the tensions and individual relation matrices should be stored in the matrix and dumped at the end computation (required for submatrix recomputation)
        :param sourced: if true, all the raltions will be looked up and not computed. Useful for the retrieval of subcirculation group, but requires the UP2voltage_and_circ to be pre-filled
        :param incremental: if True, all the circulation computation will be added to the existing ones. Usefull for the computation of particularly big systems with intermediate dumps
        :param cancellation: divides the final current by #Nodes**2/2, i.e. makes the currents comparable between circulation systems of differnet  sizes.
        :param sparse_samples: if set to an integer the sampling will be sparse and not dense, i.e. instead of compution
                                for each node pair, only an estimation will be made, equal to coputing sparse_samples association with other randomly chosen nodes
        :type sparse_samples: int
        :return: adjusted conduction system
        """
        if not incremental or self.current_accumulator == np.zeros((2,2)):
            self.current_accumulator = lil_matrix(self.Conductance_Matrix.shape)
            self.UP2UP_voltages = {}
            self.node_current = defaultdict(float)
            if not sourced:
                self.UP2voltage_and_circ = {}

        if sparse_samples:
            current_accumulator = CR.sample_pairwise_flow(self.Conductance_Matrix,
                                                                                   [self.NodeID2MatrixNumber[UP] for UP in self.analytic_Uniprots],
                                                                                   resamples=sparse_samples,
                                                                                   cancellation=cancellation)
        else:
            current_accumulator, UP_pair2voltage_current = CR.get_better_pairwise_flow(self.Conductance_Matrix,
                                                                                   [self.NodeID2MatrixNumber[UP] for UP in self.analytic_Uniprots],
                                                                                   cancellation=cancellation,
                                                                                   memoized=memoized,
                                                                                   memory_source=self.UP2voltage_and_circ)
            # self.UP2voltage_and_circ.update(UP_pair2voltage_current)
            self.UP2UP_voltages.update(dict((key, val1) for key, (val1, val2) in UP_pair2voltage_current.iteritems()))

        if incremental:
            self.current_accumulator = self.current_accumulator + current_accumulator
        else:
            self.current_accumulator = current_accumulator

        if memoized:
            self.dump_memoized()

        index_current = CR.get_current_through_nodes(self.current_accumulator)
        print current_accumulator.shape
        self.node_current.update(dict((self.MatrixNumber2NodeID[idx], val) for idx, val in enumerate(index_current)))


    def format_Node_props(self, node_current, limit = 0.01):
        """
        Formats the nodes for the analysis by in the knowledge_access_analysis module

        :param node_current: Current through the entity
        :param limit: hard limit to filter out the GO terms with too little current (co   mpensates the minor currents in the gird)
        :return: {Entity:[node current, node degree]}
        """
        charDict = {}
        limcurr = max(node_current.values())*limit
        for NodeID, i in self.NodeID2MatrixNumber.iteritems():
            if node_current[NodeID] > limcurr:
                charDict[NodeID] = [ node_current[NodeID],
                                     self.Conductance_Matrix[i, i]]
        return charDict


    def export_conduction_system(self):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format and flushes it into a file that
        can be viewed with Gephi

        :raise Warning:
        """

        if self.incomplete_compute:
            raise Warning('Computation of the information circulation was not complete, most likely due to the sampling')

        nodecharnames = ['Current', 'Type', 'Legacy_ID', 'Names', 'Degree', 'Source']
        nodechartypes = ['DOUBLE', 'VARCHAR', 'VARCHAR', 'VARCHAR', 'DOUBLE', 'DOUBLE']
        charDict = {}

        print self.Conductance_Matrix.shape
        print max(self.MatrixNumber2NodeID.iterkeys())

        for NodeID in self.node_current.iterkeys():
            i = self.NodeID2MatrixNumber[NodeID]
            charDict[NodeID] = [ str(self.node_current[NodeID]),
                             self.ID2Type[NodeID],
                             self.ID2LegacyId[NodeID],
                             self.ID2displayName[NodeID].replace(',', '-'),
                             str(self.Conductance_Matrix[i, i]),
                             str(float(int(NodeID in self.analytic_Uniprots)))]

        GDF_exporter = GDF_export_Interface(target_fname = Outputs.Interactome_GDF_output, field_names = nodecharnames,
                                            field_types = nodechartypes, node_properties_dict = charDict,
                                            mincurrent = 0.01, Idx2Label = self.MatrixNumber2NodeID,
                                            Label2Idx = self.NodeID2MatrixNumber,
                                            current_Matrix = self.current_accumulator)
        GDF_exporter.write()


    def export_subsystem(self, UP_system, UP_subsystem):
        """
        Exports the subsystem of Uniprots and circulation between them based on a larger precalculated system.This is
        possible only of the memoization parameter was on during the execution of "build_extended_circulation_system()"
        function execution.

        :param UP_system: The set of uniprots for which the larger system was calculated
        :param UP_subsystem: the set of Uniprots we are interested in
        :raise Exception: if the set of uniprots for which the larger system was calculated doesn't correspond to what
                            is stored in the dumps
        """
        current_recombinator = self.undump_memoized()
        if not set(UP_system) == set(pickle.loads(current_recombinator['UPs'])):
            raise Exception('Wrong UP system re-analyzed')
        self.UP2voltage_and_circ = pickle.loads(current_recombinator['voltages'])
        self.set_Uniprot_source(UP_subsystem)
        self.build_extended_conduction_system(memoized = False, sourced = True)
        self.export_conduction_system()


    def randomly_sample(self, samples_size, samples_each_size, sparse_rounds=False, chromosome_specific=False, memoized=False, No_add=False):
        """
        Randomly samples the set of Uniprots used to create the model. This is the null model creation routine


        :param samples_size: list of numbers of uniprots we would like to create the model for
        :param samples_each_size: how many times we would like to sample each unirot number
        :param sparse_rounds:  if we want to use sparse sampling (usefull in case of large uniprot sets),
        we would use this option
        :type sparse_rounds: int
        :param chromosome_specific: if we want the sampling to be chromosome-specific, set this parameter to the
        number of chromosome to sample from
        :type chromosome_specific: int
        :param memoized: if set to True, the sampling would be rememberd for export. Usefull in case of the chromosome comparison
        :param No_add: if set to True, the result of sampling will not be added to the database of samples. Usefull if re-running tests with similar parameters several times.
        :raise Exception: if the number of items in the samples size ann saples_each size are different
        """
        if not len(samples_size) == len(samples_each_size):
            raise Exception('Not the same list sizes!')

        self_connectable_UPs = [NodeID for NodeID, idx in self.NodeID2MatrixNumber.iteritems() if idx<(self.Conductance_Matrix.shape[0]-1)]

        if self.background is not None:
            self_connectable_UPs = list(set(self_connectable_UPs).intersection(set(self.background)))

        if chromosome_specific:
            self_connectable_UPs = list(set(self_connectable_UPs).intersection(set(self.Chrom2UP[str(chromosome_specific)])))

        for sample_size, iterations in zip(samples_size, samples_each_size):
            for i in range(0, iterations):
                shuffle(self_connectable_UPs)
                analytics_UP_list = self_connectable_UPs[:sample_size]
                self.set_Uniprot_source(analytics_UP_list)
                self.build_extended_conduction_system(memoized=memoized, sourced=False, sparse_samples=sparse_rounds)

                md5 = hashlib.md5(json.dumps( sorted(analytics_UP_list), sort_keys=True)).hexdigest()

                # print 'debug1 \t', self.UP2UP_voltages

                if not No_add:
                    Interactome_rand_samp.insert({'UP_hash' : md5,
                                     'sys_hash' : self._MD5hash(),
                                     'size' : sample_size,
                                     'chrom': str(chromosome_specific),
                                     'sparse_rounds': sparse_rounds,
                                     'UPs' : pickle.dumps(analytics_UP_list),
                                     'currents' : pickle.dumps((self.current_accumulator, self.node_current)),
                                     'voltages' : pickle.dumps(self.UP2UP_voltages)})

                print 'Random ID: %s \t Sample size: %s \t iteration: %s\t compops: %s \t time: %s ' %(self.r_ID,
                        sample_size, i, "{0:.2f}".format(sample_size*(sample_size-1)/2/self._time()), self.pretty_time())


    def map_UPs_to_chromosomes(self):
        """
        Maps the Uniprot Ids to chromosomes of the organism and reversely.

        """
        resdict = {}
        # TODO: add the possibility that Chromosome_source wasn't properly declared or is a mapping to a void declaration
        #          if this is the case,
        for fle in os.listdir(Chromosome_source):
            if Chromosome_file_filter in fle:
                source_fle = Chromosome_source+'/'+fle
                id = fle.split('.')[0].split(Chromosome_file_filter)[1]
                resdict[id]=open(source_fle).read()

        for i, (key, val) in enumerate(Mat_gter.ID2LegacyId.iteritems()):
            for id, textblock in resdict.iteritems():
                if val in textblock:
                    self.UP2Chrom[key] = id
                    self.Chrom2UP[id].append(key)

        self.Chrom2UP = dict(self.Chrom2UP)



if __name__ == "__main__":
    Mat_gter = MatrixGetter(True, True)
    Mat_gter.full_rebuild ()
    # Mat_gter.fast_load()
    print Mat_gter.pretty_time()

    # test_set = ['147875', '130437', '186024', '100154', '140777', '100951', '107645', '154772']
    #
    # test2 = ['55618', '55619', '55616', '55614', '55615', '55 612', '55613', '55342', '177791', '126879']
    #
    # Mat_gter.set_Uniprot_source(test_set)
    # Mat_gter.export_subsystem(test_set, test2)
    # Mat_gter.build_extended_conduction_system()
    # Mat_gter.export_conduction_system()

    # Mat_gter.randomly_sample([100,250],[5,5], sparse_rounds=10)


