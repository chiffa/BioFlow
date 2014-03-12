__author__ = 'ank'
"""
This module contains all the routines that are respojnsible for pulling
the matrixes out of the neo4j graph and processing them.
"""
import itertools
import numpy as np
import pickle
from copy import copy
from time import time

from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from PolyPharma.configs import edge_type_filters, Adjacency_Martix_Dict, Conductance_Matrix_Dict, Dumps
from PolyPharma.neo4j_analyzer.IO_Routines import reaction_participant_getter, expand_from_seed, Erase_custom_fields

# TODO: change the behavior of HiNT propagation to a several stages propagation
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

    # # List of all the reaction types present in the DatabaseGraph that will be
    # # used as roots to build the interaction network (not all nodes are necessary
    # # within the connex part of the graph)
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
        self.ID2Type = {}
        self.ID2Localization = {}
        self.Uniprots = []
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


    def write_to_csv(self, filename, array):
        """
        Dumps a numpy array to csv

        :param filename: location of dumping
        :type filename: str
        :param array: array to dump
        :type array: numpy.array
        """
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
        """
        Undumps a pickled object

        :param dump_filename: filename from which undump
        :type dump_filename: str
        :return: the undumped object
        :rtype: object
        """
        DF = file(dump_filename, 'r')
        return pickle.load(DF)


    def dump_Matrices(self):
        """
        dumps self.Ajacency_Matrix and self.Conductance_Matrix
        """
        self.dump_object(Dumps.ValMat, self.Ajacency_Matrix)
        self.dump_object(Dumps.ConMat, self.Conductance_Matrix)


    def undump_Matrices(self):
        """
        undumps self.Ajacency_Matrix and self.Conductance_Matrix
        """
        self.Ajacency_Matrix = self.undump_object(Dumps.ValMat)
        self.Conductance_Matrix = self.undump_object(Dumps.ConMat)


    def dump_Eigens(self):
        """
        dumps self.adj_eigenvals and self.Conductance_Matrix and writes them to csv
        """
        self.write_to_csv(Dumps.eigen_VaMat, self.adj_eigenvals)
        self.write_to_csv(Dumps.eigen_ConMat, self.cond_eigenvals)
        self.dump_object(Dumps.val_eigen, (self.adj_eigenvals, self.adj_eigenvects))
        self.dump_object(Dumps.cond_eigen, (self.cond_eigenvals, self.cond_eigenvects))


    def undump_Eigens(self):
        """
        undumps self.adj_eigenvals and self.Conductance_Matrix
        """
        self.adj_eigenvals, self.adj_eigenvects = self.undump_object(Dumps.val_eigen)
        self.cond_eigenvals, self.cond_eigenvects = self.undump_object(Dumps.cond_eigen)



    def dump_Maps(self):
        """
        dumps all the elements required for the mapping between the types and ids of database entries and matrix columns
        """
        self.dump_object(Dumps.matrix_corrs,
                         (self.NodeID2MatrixNumber, self.MatrixNumber2NodeID,
                          self.ID2displayName, self.ID2Type, self.ID2Localization,
                          self.Uniprots, self.Uniprot_attachments, self.Uniprot_Mat_idxs))


    def undump_Maps(self):
        """
        undumps all the elements required for the mapping between the types and ids of database entries and matrix columns
        """
        self.NodeID2MatrixNumber, self.MatrixNumber2NodeID, self.ID2displayName, self.ID2Type,\
        self.ID2Localization, self.Uniprots, self.Uniprot_attachments, self.Uniprot_Mat_idxs = self.undump_object(Dumps.matrix_corrs)

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

        # self.dump_object('fixture.dump',(self.ReactLinks,self.InitSet,self.GroupLinks, self.GroupSet, self.SecLinks, self.SecSet, self.UP_Links, self.UPSet, self.HiNT_Links, self.FullSet, self.Super_Links, self.ExpSet))


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

        for ID in self.Highest_Set:
            self.NodeID2MatrixNumber[ID] = counter
            self.MatrixNumber2NodeID[counter] = ID
            Vertex = DatabaseGraph.vertices.get(ID)
            self.ID2displayName[ID] = Vertex.displayName
            self.ID2Type[ID] = Vertex.element_type
            if Vertex.element_type == "UNIPROT":
                self.Uniprots.append(ID)
                self.Uniprot_Mat_idxs.append(counter)
            if Vertex.localization is not None:
                self.ID2Localization[ID] = request_location(LocationBufferDict, Vertex.localization)
            counter += 1


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
        print "Marking of %s nodes for connexity was done in %s" % (str(ln), str(self.time()))


    def full_rebuild(self):
        """
        Performs the initial loading routines that set up in place the sytem of dump files to allow fast loading in the
        future
        """

        CA  = self.Connexity_Aware
        self.Connexity_Aware = False

        self.create_val_matrix()

        # TODO: hard connexity reset
        # TODO: better selection of UP so that they are adherent to Reactome_nodes
        Erase_custom_fields()
        self.Write_Connexity_Infos()

        self.compute_Uniprot_Attachments()
        self.Connexity_Aware = CA
        self.create_val_matrix()
        self.get_eigenspectrum(100)

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
        self.undump_Maps()
        self.Uniprot_Mat_idxs = []
        for SP_Id in self.Uniprots:
            Node=DatabaseGraph.UNIPORT.get(SP_Id)
            if Node.main_connex:
                self.Uniprot_Mat_idxs.append(self.NodeID2MatrixNumber[SP_Id])
        print len(self.Uniprot_Mat_idxs)
        self.dump_Maps()


if __name__ == "__main__":
    Mat_gter = MatrixGetter(True, True)
    # Mat_gter.hacky_corr()
    Mat_gter.full_rebuild ()
    # Mat_gter.fast_load()
    Mat_gter.get_eigenspectrum(100)
    Mat_gter.dump_Eigens()