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
from PolyPharma.configs import IDFilter
from copy import copy
from time import time
import itertools
import numpy as np
# TODO: remove the sadistic line below !!!!!!!!!!!!!!!!!!!!
import json
# TODO: remove the sadistic line above !!!!!!!!!!!!!!!!!!!!
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import connected_components

# TODO: switch from filtering to annotating as part of Cross-ref set


# Refers to the groups of links between the nodes that should be treated in the same manner
edge_type_filter0_1 = ["is_part_of_collection"]                    # Group relation group
edge_type_filter0_2 = ["is_same"]                                  # Same relation group
edge_type_filter1_1 = ["is_Catalysant", "is_reaction_particpant"]  # Reaction relation group
edge_type_filter1_2 = ["is_part_of_complex", "is_Regulant"]        # Contact_interaction relation group
edge_type_filter1_3 = ["is_interacting"]                           # Contact_interaction relation group
edge_type_filter2 = ["is_possibly_same"]                           # possibly_same relation group

#TODO: move the coefficients to the configs file

# Coefficients values for the value_Matrix
DfactorDict = {"Group":0.5,
             "Same":1,
             "Reaction":0.33,
             "Contact_interaction":0.33,
             "possibly_same":0.1,
             }

# Coefficients values for the conductance_Matrix
ConductanceDict = {"Group":0.5,
             "Same":100,
             "Reaction":1,
             "Contact_interaction":1,
             "possibly_same":0.1,
             }

# List of all the reaction types present in the DatabaseGraph that will be
# used as roots to build the interaction network (not all nodes are necessary
# within the connex part of the graph)
ReactionsList = [DatabaseGraph.TemplateReaction, DatabaseGraph.Degradation, DatabaseGraph.BiochemicalReaction]

class Dumps(object):
    matrix_LS = 'dump5.dump'
    matrix_corrs = 'dump2.dump'
    eigen_VaMat = 'eigen_valmat.csv'
    eigen_ConMat = 'eigen_conmat.csv'
    val_eigen = 'pickleDump.dump'
    cond_eigen = 'pickleDump_0_5.dump'
    ValMat = 'pickleDump3.dump'
    ConMat = 'pickleDump4.dump'
    UniP_att = 'UP_Attach.dump'


def getMatrix(decreaseFactorDict, numberEigvals, FastLoad, ConnexityAwareness, full_impact):

    def build_correspondances(IDSet, Rapid):
        """ Maps Node IDs to matrix row/column indexes; """

        def request_location(LocationBufferDict, location):
            """Buffered lookup of location"""
            location = str(location)
            if location in LocationBufferDict.keys():
                return LocationBufferDict[location]
            else:
                generator = DatabaseGraph.Location.index.lookup(ID = location)
                if generator != None:
                    for elt in generator:
                        LocationBufferDict[location] = str(elt.displayName)
                        return str(elt.displayName)

        NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, Uniprots, ID2Localization = ({},{},{},{},[],{})
        counter = 0
        LocationBufferDict = {}

        if Rapid:
            DF = file(Dumps.matrix_corrs, 'r')
            NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = json.load(DF)
            return NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots

        for ID in IDSet:
            NodeID2MatrixNumber[ID] = counter
            MatrixNumber2NodeID[counter] = ID
            Vertex = DatabaseGraph.vertices.get(ID)
            ID2displayName[ID] = Vertex.displayName
            ID2Type[ID] = Vertex.element_type
            if Vertex.element_type == "UNIPROT":
                Uniprots.append(ID)
            if Vertex.localization != None:
                ID2Localization[ID] = request_location(LocationBufferDict, Vertex.localization)
            counter += 1
        DF = file(Dumps.matrix_corrs,'w')
        json.dump((NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots), DF)
        DF.close()
        return NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots

    # noinspection PyTypeChecker
    def Full_load(Connexity_Aware, init_time):

        def get_Reaction_blocks(Connexity_Aware):
            '''
            Recovers the blocks if interaction that are due to a common set of reactions
            for the elements. They will be used as roots to build the complete interaction
            tree later on.

            :param Connexity_Aware: if this parameter is set for true, only the elemets that are within the major connexity
                                    graph will be loaded
            :type Connexity_Aware: boolean

            .. warning::
            Do not set Connexity_Aware to True on the first run or before performing the Write_Connexity_Infos routine

            :returns: list of lists of NodeIDs -- Reagent_Clusters: Clusters of reagents that interact together through
                            a common reaction
            :returns: List of NodeIDs -- Seeds: NodeIDs of the nodes that were encountered, used to perform a further
                            expansion regarding the other links (Complex participation, contact, etc...)
            :returns: int -- count: number of interaction clusters
            '''
            ReagentClusters = []
            Seeds = set()
            count = 0
            for ReactionType in ReactionsList:
                for Reaction in ReactionType.get_all():
                    if Reaction == None:
                        continue
                    LocalList = []
                    for edge_type in edge_type_filter1_1:
                        if Reaction.bothV(edge_type) == None:
                            continue
                        for elt in Reaction.bothV(edge_type):
                            Connex = True

                            if Connexity_Aware:
                                Connex = False
                                if  elt.custom != None and "Main_Connex" in elt.custom:
                                    Connex = True

                            ID = str(elt).split('/')[-1][:-1]
                            if ID not in IDFilter and Connex:
                                LocalList.append(ID)
                                count += 1

                    if len(LocalList) > 1:
                        ReagentClusters.append(copy(LocalList))
                        Seeds.update(LocalList)

            return ReagentClusters, Seeds, count

        def get_expansion(SubSeed, edge_type_filter):
            '''
            Recovers all the nodes reached from the SubSeed according to a the relations listed
            in edge_type_filter

            :param SubSeed: List of NodeIDs that serve as a root for searching further relations
            :type SubSeed: List of NodeIDs

            :param edge_type_filter: type of relations according to which the graph will be explored
            :type edge_type_filter: relations of the type bulbsGraph.Relationtype (these are well python objects)

            :return Clusters: Dictionary of lists, where key is the NodeID from the SubSeed
            :return SuperSet: List of NodeIDs attained from the SubSet by the relation types from the edge_type_filter
            '''
            Clusters = {}
            SuperSeed = set()
            SuperSeed.update(SubSeed)
            count = 0
            # noinspection PyTypeChecker
            for element in SubSeed:
                SeedNode = DatabaseGraph.vertices.get(element)
                LocalList = []
                for edge_type in edge_type_filter:
                    if SeedNode.bothV(edge_type) != None:
                        for elt in SeedNode.bothV(edge_type):
                            ID = str(elt).split('/')[-1][:-1]
                            if ID not in IDFilter:
                                LocalList.append(ID)
                                SuperSeed.add(ID)
                                count += 1
                if len(LocalList) > 0:
                    Clusters[element] = copy(LocalList)
            return Clusters, SuperSeed, count

        def characterise(Links, Group, c, init_time, t_time):
            print len(Links), len(Group), c
            print time() - init_time, time() - t_time
            return time()

        t = time()

        ReactLinks, InitSet, c = get_Reaction_blocks(ConnexityAwareness)
        t = characterise(ReactLinks, InitSet, c, init_time, t)

        GroupLinks, GroupSet, c = get_expansion(InitSet, edge_type_filter0_1)
        t = characterise(GroupLinks, GroupSet, c, init_time, t)

        SecLinks, SecSet, c = get_expansion(GroupSet, edge_type_filter1_2)
        t = characterise(SecLinks, SecSet, c, init_time, t)

        for i in range(0,5):
            SecLinks2, SecSet2, c = get_expansion(SecSet, edge_type_filter1_2)
            t = characterise(SecLinks2, SecSet2, c, init_time, t)
            SecSet = SecSet2
            SecLinks = SecLinks2

        UP_Links, UPSet, c = get_expansion(SecSet, edge_type_filter0_2)
        t = characterise(UP_Links, UPSet, c, init_time, t)

        HiNT_Links, FullSet, c = get_expansion(UPSet, edge_type_filter1_3)
        t = characterise(HiNT_Links, FullSet, c, init_time, t)

        Super_Links, ExpSet, c = get_expansion(FullSet, edge_type_filter2)
        characterise(Super_Links, ExpSet, c, init_time, t)

        return ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet

    def fast_insert(element, index_type):
        # TODO: improve to pump from two different dictionnaries for ValueMatrix and ConductanceMatrix.
        ValueMatrix[element[0], element[1]] = min(ValueMatrix[element[0], element[1]] + decreaseFactorDict[index_type], 1)
        ValueMatrix[element[1], element[0]] = min(ValueMatrix[element[1], element[0]] + decreaseFactorDict[index_type], 1)
        ConductanceMatrix[element[0], element[1]] = ConductanceMatrix[element[0], element[1]] - decreaseFactorDict[index_type]
        ConductanceMatrix[element[1], element[0]] = ConductanceMatrix[element[1], element[0]] - decreaseFactorDict[index_type]
        ConductanceMatrix[element[1], element[1]] = ConductanceMatrix[element[1], element[1]] + decreaseFactorDict[index_type]
        ConductanceMatrix[element[0], element[0]] = ConductanceMatrix[element[0], element[0]] + decreaseFactorDict[index_type]


    init_time = time()

    ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks,\
    SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet = ({},[],{},[],{},[],{},[],{},[],{},[])

    if not FastLoad:
        ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet = Full_load(ConnexityAwareness, init_time)
        DF = file(Dumps.matrix_LS, 'w') # previously 'dump5.dump', 'w'
        json.dump((ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet),DF)
        DF.close()

    else:
        DF = file(Dumps.matrix_LS, 'r')
        ReactLinks, InitSet, GroupLinks, GroupSet, SecLinks, SecSet, UP_Links, UPSet, HiNT_Links, FullSet, Super_Links, ExpSet = json.load(DF)

    Highest_Set = FullSet

    if full_impact:
        Highest_Set = ExpSet

    print "building correspondances", time() - init_time

    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = build_correspondances(Highest_Set, FastLoad)

    print "building the ValMatrix", time() - init_time

    loadLen = len(Highest_Set)
    ValueMatrix = lil_matrix((loadLen, loadLen))
    ConductanceMatrix = lil_matrix((loadLen, loadLen))

    # the difference between ReactLinks and all the others are that React_Links are eqally conected groups,
    # whereath all the other are made by performing inclusion into an existing element and are of form dict[root]=[leaves]
    # And thus we can actually perform the correlation isertion through the same function, as long as we iterate on non-
    # repeated pairs.

    t = time()

    for group in ReactLinks:
        for elt in itertools.combinations(group, 2):
            element = (NodeID2MatrixNumber[elt[0]], NodeID2MatrixNumber[elt[1]])
            fast_insert(element, "Reaction")

    for key in GroupLinks.keys():
        for val in GroupLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            fast_insert(element, "Group")

    for key in SecLinks.keys():
        for val in SecLinks[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            fast_insert(element, "Contact_interaction")

    for key in UP_Links.keys():
        for val in UP_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            fast_insert(element, "Same")

    for key in HiNT_Links.keys():
        for val in HiNT_Links[key]:
            element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
            fast_insert(element, "Contact_interaction")

    if full_impact:
        for key in Super_Links.keys():
            for val in Super_Links[key]:
                element=(NodeID2MatrixNumber[key],NodeID2MatrixNumber[val])
                fast_insert(element, "possibly_same")


    print "entering eigenvect computation", time()-init_time, time()-t

    eigenvals, eigenvects = eigsh(ValueMatrix, numberEigvals)
    eigenvals2, eigenvects2 = eigsh(ConductanceMatrix, numberEigvals)

    print eigenvals
    print '<======================>'
    print eigenvals2
    print '<======================>'
    print np.all(eigsh(ConductanceMatrix)[0] > 0)

    DF = file(Dumps.eigen_VaMat, 'w')
    DF.write(eigenvals)
    DF.close()

    DF = file(Dumps.eigen_ConMat, 'w')
    DF.write(eigenvals2)
    DF.close()

    DF = file(Dumps.val_eigen, 'w')
    json.dump((eigenvals, eigenvects), DF)
    DF.close()

    DF = file(Dumps.cond_eigen, 'w')
    json.dump((eigenvals2,eigenvects2), DF)
    DF.close()

    DF = file(Dumps.ValMat, 'w')
    json.dump(ValueMatrix, DF)
    DF.close()

    DF=file(Dumps.ConMat, 'w')
    json.dump(ConductanceMatrix, DF)
    DF.close()

    print time()-init_time, time()-t


def Write_Connexity_Infos():
    """
        Writes the infos about connexity of different components of a graph. This execution is the main
        reason for the existance of the "Value" Matrix.
    """
    #TODO: build unittest to see if it works correctly.
    ValueMatrix = json.load(file(Dumps.ValMat, 'r'))
    CCOmps = connected_components(ValueMatrix, directed = False)
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = json.load(file(Dumps.matrix_corrs,'r'))

    counters = np.zeros((CCOmps[0], 1))
    for elt in range(0, len(CCOmps[1])):
        counters[CCOmps[1][elt], 0] += 1
    major_Index = np.argmax(counters)

    ln = len(CCOmps[1])
    for i in range(0, len(CCOmps[1])):
        print 'running,', i, ln
        if CCOmps[1][i] == major_Index:
            Node = DatabaseGraph.vertices.get(MatrixNumber2NodeID[i])
            Node.custom = 'Main_Connex'
            Node.save()

def Erase_Additional_Infos():
    """
        Resets the .costum field of all the Nodes on which we have iterated here. Usefull to perform
        After node set or node connectivity were modfied.
    """
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = json.load(file(Dumps.matrix_corrs,'r'))

    for NodeID in NodeID2MatrixNumber.keys():
        Node = DatabaseGraph.vertices.get(NodeID)
        Node.custom = ''
        Node.save()

def compute_Uniprot_Attachments():
    """
        Attaches the Uniprots to the proteins from the reactome, allowing to combine the information circulation values
    """
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = json.load(file(Dumps.matrix_corrs,'r'))
    Uniprot_Attach = {}
    #TODO: Build a test to see if this one is still truly required

    for SP_Node_ID in Uniprots:
        Vertex = DatabaseGraph.UNIPORT.get(SP_Node_ID)
        Generator = Vertex.bothV("is_same")
        if Generator != None:
            Uniprot_Attach[SP_Node_ID] = []
            for item in Generator:
                ID = str(item).split('/')[-1][:-1]
                Uniprot_Attach[SP_Node_ID].append(ID)
        print SP_Node_ID, 'processed', len(Uniprot_Attach[SP_Node_ID])
    json.dump(Uniprot_Attach,file(Dumps.UniP_att,'w'))
    return Uniprot_Attach


def Perform_Loading_Routines():
    '''
    Performs the initial loading routines that set up in place the sytem of dump files and co
    for the main algorithms to kick in seamlessly
    '''

    getMatrix(DfactorDict, 1, False, False, False)
    Erase_Additional_Infos()
    Write_Connexity_Infos()
    getMatrix(DfactorDict, 100, False, True, False)
    compute_Uniprot_Attachments()

if __name__ == "__main__":
    pass