__author__ = 'ank'
"""
This module contains all the routines that are respojnsible for pulling
the matrixes out of the neo4j graph and processing them

The general idea is to build up
- a value matrix that references only the connexions between distinct nodes
- a conductance matrix that refers the conductance between the nodes and the
self-referenced conductances
"""

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from PolyPharma.configs import edge_type_filters, Adjacency_Martix_Dict, Conductance_Matrix_Dict, Dumps
from copy import copy
from time import time
import itertools
import numpy as np
import pickle
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components
from IO_Routines import reaction_participant_getter,expand_from_seed



class MatrixGetter(object):


    # # List of all the reaction types present in the DatabaseGraph that will be
    # # used as roots to build the interaction network (not all nodes are necessary
    # # within the connex part of the graph)
    ReactionsList = [DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]


    def __init__(self, FastLoad, Connexity_Aware, full_impact):
        """
        .. code-block:: python
        >>> data=['some cool stuff']
        """

        self.FastLoad = FastLoad
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
        self.ID2Type = {}
        self.ID2Localization ={}
        self.Uniprots = []

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

        self.main_connexity_set_IDs = []

        self.Matrix_reality = Dumps.matrix_corrs


    def write_to_csv(self, filename, array):
        DF = file(filename, 'w')
        DF.write(array)
        DF.close()


    def dump_object(self, dump_filename, object_to_dump):
        """
        Shorcut for pickling & dumping behavior

        :param dump_filename: filename where the object will be dumped
        :type dump_filename: str

        :param object_to_dump: object to be pickled and dumped
        :type object_to_dump: pickable object
        """
        DF = file(dump_filename, 'w')
        pickle.dump(object_to_dump, DF)
        DF.close()


    def undump_object(self, dump_filename):
        DF = file(dump_filename, 'r')
        return pickle.load(DF)


    def dump_Matrices(self):
        self.dump_object(Dumps.ValMat, self.Ajacency_Matrix)
        self.dump_object(Dumps.ConMat, self.Conductance_Matrix)


    def undump_Matrices(self):
        self.ValueMatrix = self.undump_object(Dumps.ValMat)
        self.Conductance_Matrix = self.undump_object(Dumps.ConMat)


    def dump_Eigens(self):
        self.write_to_csv(Dumps.eigen_VaMat, self.adj_eigenvals)
        self.write_to_csv(Dumps.eigen_ConMat, self.cond_eigenvals)
        self.dump_object(Dumps.val_eigen, (self.adj_eigenvals, self.adj_eigenvects))
        self.dump_object(Dumps.cond_eigen, (self.cond_eigenvals, self.cond_eigenvects))


    def undump_Eigens(self):
        self.adj_eigenvals, self.adj_eigenvects = self.undump_object(Dumps.val_eigen)
        self.cond_eigenvals, self.cond_eigenvects = self.undump_object(Dumps.cond_eigen)



    def dump_Maps(self):
        self.dump_object(Dumps.matrix_corrs,
                         (self.NodeID2MatrixNumber, self.MatrixNumber2NodeID,
                          self.ID2displayName, self.ID2Type, self.ID2Localization, self.Uniprots))


    def undump_Maps(self):
        self.NodeID2MatrixNumber, self.MatrixNumber2NodeID, self.ID2displayName, self.ID2Type, self.ID2Localization, self.Uniprots = self.undump_object(Dumps.matrix_corrs)


    def dump_Main_connex_set(self):
        self.dump_object(Dumps.Main_Connex_group,self.main_connexity_set_IDs)


    def undump_Main_connex_set(self):
        self.main_connexity_set_IDs = self.undump_object(Dumps.Main_Connex_group)


    def time(self):
        it, pt = (round(time() - self.init_time), round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (int(it) / 60, it % 60, int(pt) / 60, pt % 60)
        self.partial_time = time()
        return pload


    def full_load_LS(self):

        def get_Reaction_blocks():
            '''
            Recovers the blocks if interaction that are due to a common set of reactions
            for the elements. They will be used as roots to build the complete interaction
            tree later on.

            '''
            ReagentClusters = []
            Seeds = set()
            count = 0
            for ReactionType in self.ReactionsList:
                for Reaction in ReactionType.get_all():
                    if Reaction == None:
                        continue
                    LocalList, cnt = reaction_participant_getter(Reaction, self.Connexity_Aware)
                    count += cnt

                    if len(LocalList) > 1:
                        ReagentClusters.append(copy(LocalList))
                        Seeds.update(LocalList)

            return ReagentClusters, Seeds, count

        def get_expansion(SubSeed, edge_type_filter):
            '''
            Recovers all the nodes reached from the SubSeed according to a the relations listed
            in edge_type_filter

            :param SubSeed: set of NodeIDs that serve as a root for searching further relations
            :type SubSeed: set of NodeIDs

            :param edge_type_filter: type of relations according to which the graph will be explored
            :type edge_type_filter: relations of the type bulbsGraph.Relationtype (these are well python objects)

            :return Clusters: Dictionary of lists, where key is the NodeID from the SubSeed
            :return SuperSet: List of NodeIDs attained from the SubSet by the relation types from the edge_type_filter
            '''
            Clusters = {}
            SuperSeed = set()
            SuperSeed.update(SubSeed)
            count = 0
            for element in SubSeed:
                LocalList, count_increase = expand_from_seed(element, edge_type_filter)
                if len(LocalList) > 0:
                    Clusters[element] = copy(LocalList)
                    SuperSeed.update(LocalList)
                    count+=count_increase
            return Clusters, SuperSeed, count


        def characterise(name, Links, Group, c):
            print '===========>', name, '<==========='
            print len(Links), len(Group), c
            print self.time()

        ################################################################################################################

        self.ReactLinks, self.InitSet, c = get_Reaction_blocks()
        characterise('Reactions', self.ReactLinks, self.InitSet, c)

        self.GroupLinks, self.GroupSet, c = get_expansion(self.InitSet, edge_type_filters["Group"])
        characterise('Groups', self.GroupLinks, self.GroupSet, c)

        self.SecLinks, self.SecSet, c = get_expansion(self.GroupSet, edge_type_filters["Contact_interaction"])
        characterise('Secondary Links', self.SecLinks, self.SecSet, c)

        for i in range(0,5):
            SecLinks2, SecSet2, c = get_expansion(self.SecSet, edge_type_filters["Contact_interaction"])
            self.SecSet = SecSet2
            self.SecLinks = SecLinks2
            characterise('Secondary Links '+str(i)+' ', self.SecLinks, self.SecSet, c)

        self.UP_Links, self.UPSet, c = get_expansion(self.SecSet, edge_type_filters["Same"])
        characterise('Uniprot Links', self.UP_Links, self.UPSet, c)

        self.HiNT_Links, self.FullSet, c = get_expansion(self.UPSet, edge_type_filters["HiNT_Contact_interaction"])
        characterise('HiNT Links', self.HiNT_Links, self.FullSet, c)

        self.Super_Links, self.ExpSet, c = get_expansion(self.FullSet, edge_type_filters["possibly_same"])
        characterise('Looks_similar Links', self.Super_Links, self.ExpSet, c)

        self.dump_object('fixture.dump',(self.ReactLinks,self.InitSet,self.GroupLinks, self.GroupSet, self.SecLinks, self.SecSet, self.UP_Links, self.UPSet, self.HiNT_Links, self.FullSet, self.Super_Links, self.ExpSet))


    def map_rows_to_names(self,):
        """ Maps Node IDs to matrix row/column indexes; """

        def request_location(LocationBufferDict, location):
            """Just a Buffered lookup of location"""
            location = str(location)
            if location in LocationBufferDict.keys():
                return LocationBufferDict[location]
            else:
                generator = DatabaseGraph.Location.index.lookup(ID = location)
                if generator != None:
                    for elt in generator:
                        LocationBufferDict[location] = str(elt.displayName)
                        return str(elt.displayName)

        ######################################################################################################################


        counter = 0
        LocationBufferDict = {}

        for ID in self.Highest_Set:
            self.NodeID2MatrixNumber[ID] = counter
            self.MatrixNumber2NodeID[counter] = ID
            Vertex = DatabaseGraph.vertices.get(ID)
            self.ID2displayName[ID] = Vertex.displayName
            self.ID2Type[ID] = Vertex.element_type
            if Vertex.element_type == "UNIPROT":
                self.Uniprots.append(ID)
            if Vertex.localization != None:
                self.ID2Localization[ID] = request_location(LocationBufferDict, Vertex.localization)
            counter += 1


    def fast_row_insert(self, element, index_type):
        """
        performs an correct insertion of an edge to the matrix.
        """

        self.Ajacency_Matrix[element[0], element[1]] = min(self.Ajacency_Matrix[element[0], element[1]] + Adjacency_Martix_Dict[index_type], 1)
        self.Ajacency_Matrix[element[1], element[0]] = min(self.Ajacency_Matrix[element[1], element[0]] + Adjacency_Martix_Dict[index_type], 1)

        self.Conductance_Matrix[element[0], element[1]] -= Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[1], element[0]] -= Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[1], element[1]] += Conductance_Matrix_Dict[index_type]
        self.Conductance_Matrix[element[0], element[0]] += Conductance_Matrix_Dict[index_type]


    def create_val_matrix(self):
        """
        Creates the Value and Conductance matrices
        """

        self.full_load_LS()

        self.Highest_Set = self.FullSet

        if self.full_impact:
            self.Highest_Set = self.ExpSet

        print "building correspondances", self.time()

        self.map_rows_to_names()

        print "building the ValMatrix", self.time()

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

        if self.full_impact:
            for key in self.Super_Links.keys():
                for val in self.Super_Links[key]:
                    element = (self.NodeID2MatrixNumber[key], self.NodeID2MatrixNumber[val])
                    self.fast_row_insert(element, "possibly_same")


    def get_eigenspectrum(self, numberEigvals):

        if self.Conductance_Matrix == np.zeros((4,4)) or self.FullSet == []:
            raise Exception("Matrix must be pre-loaded first")

        print "entering eigenvect computation", self.time()

        self.adj_eigenvals, self.adj_eigenvects = eigsh(self.Ajacency_Matrix, numberEigvals)
        self.cond_eigenvals, self.cond_eigenvects = eigsh(self.Conductance_Matrix, numberEigvals)

        print self.adj_eigenvals
        print '<======================>'
        print self.cond_eigenvals
        print '<======================>'
        print np.all(eigsh(self.Conductance_Matrix)[0] > 0)

        print self.time()


    def get_normalized_laplacian(self, fudge = 1E-18):
        d, V = np.linalg.eigh(self.Conductance_Matrix)
        D = np.diag(1./np.sqrt(d + fudge))
        return np.dot(np.dot(V, D), V.T)


    def Write_Connexity_Infos(self):
        """
            Writes the infos about connexity of different components of a graph. This execution is the main
            reason for the existance of the "Value" Matrix.
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
                self.main_connexity_set_IDs.append(self.MatrixNumber2NodeID[i])
                Node = DatabaseGraph.vertices.get(self.MatrixNumber2NodeID[i])
                Node.custom = 'Main_Connex'
                Node.save()
        print "Marking of %s nodes for connexity was done in %s" % (str(ln), str(self.time()))


    def full_rebuild(self):
        """
        Performs the initial loading routines that set up in place the sytem of dump files and co
        """

        CA  = self.Connexity_Aware
        self.Connexity_Aware = False

        self.create_val_matrix()

        self.undump_Main_connex_set()
        self.reset_Connexity_Infos()
        self.Write_Connexity_Infos()
        self.dump_Main_connex_set()

        self.Connexity_Aware = CA
        self.create_val_matrix()
        self.get_eigenspectrum(100)


    def reset_Connexity_Infos(self):
        """
        Resets the .costum field of all the Nodes on which we have iterated here. Usefull to perform
        After node set or node connectivity were modfied.
        """

        for NodeID in self.main_connexity_set_IDs:
            Node = DatabaseGraph.vertices.get(NodeID)
            Node.save()

        print "has resetted connexity over %s nodes in %s" % (len(self.main_connexity_set_IDs), str(self.time()))


if __name__ == "__main__":
    Mat_gter = MatrixGetter(True, True, True)
    Mat_gter.full_rebuild()