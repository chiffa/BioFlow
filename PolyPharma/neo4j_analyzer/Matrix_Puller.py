__author__ = 'ank'
"""
:author: andrei

This module contains all the routines that are respojnsible for pulling 
the matrixes out of the neo4j graph and processing them
"""

if __name__ == "__main__" and __package__ is None:
    __package__ = "PolyPharma.neo4j_analyzer"

import copy
import operator
import random
import pickle
import pylab
import numpy as np
from itertools import combinations, repeat
from scipy.sparse import lil_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import diags
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky

from PolyPharma.configs import Dumps
from PolyPharma.neo4j_analyzer.Matrix_Interactome_DB_interface import MatrixGetter

ID2Aboundances = {}

# NOTE: in case of creation of a fixture to test the class functionalities, the fixture object should be imported instead of MG and fake it's objects.
MG = MatrixGetter(True, False)
MG.fast_load()


def _characterize(Objkt):
    """
    Small snipped used to debug unexpected behavior of matrix/ndarray/sparse matrix mixtures

    :param Objkt: object to characteize
    """
    print Objkt.shape, Objkt.sum(),
    print len(Objkt.nonzero()[0]),
    print type(Objkt)


def characterise_eigenvects(N_of_eigenvectors_to_characterise, Type='Adj'):
    """
    Recovers statistics over eigenvectors in order to remove excessively connected nodes.

    :param N_of_eigenvectors_to_characterise: number of eigenvectors we would like to characterise
    :param Type: if Type='Adj', the adjacency matrix eigenvectors will be processed. Otherwise a discrete laplacian will be used
    :return: Names of the elements in the 100 biggest eigenvectors in the adjacency matrix
    """
    Retname = set()
    if Type == 'Adj':
        norm = np.linalg.norm(MG.adj_eigenvects, axis = 1)
    else:
        norm = np.linalg.norm(MG.cond_eigenvects, axis = 1)
    sbigest_args = np.argsort(norm)[- N_of_eigenvectors_to_characterise:]

    for arg in reversed(sbigest_args):
        print 'Matrix row: %s, \t Row Norm: %s, \t Description: %s' % (arg,  norm[arg],  MG.get_descriptor_for_index(arg))
        Retname.add(str(MG.get_descriptor_for_index(arg)[1]))

    return Retname


def checkMatrix():
    """
    Verifyes that the adjacency matrix has no elements with unexpected values. If there are, prints them, then throws
    an exception

    :Note: this constitutes one of the possible partitioning of the matrix into functional clusters.
    :raise Exception: "Unexpected behavior among adjacency nodes!" if adjacency matrix has unexpected elements
    """
    NzeroList = MG.Ajacency_Matrix.nonzero()
    faultyList = []
    for i in range(0, len(NzeroList[0])):
        index = (int(NzeroList[0][i]), int(NzeroList[1][i]))
        value = MG.Ajacency_Matrix[index[0],index[1]]
        if value < 0.0 or value > 1.0:
            faultyList.append(index)
    print "There are %s faulty nodes out of %s" % (len(faultyList),len(NzeroList[0]))
    if len(faultyList)>0:
        print faultyList
        raise Exception("Unexpected behavior among adjacency nodes!")



def processEigenVectors(nbiggest, Type='Adj'):
    """
    Recovers "nbiggest" elements with the biggest association to a given eigenvector and returns them.

    :param nbiggest: number of biggest elements in a eigenvector to analyse
    :param Type: if Type='Adj', the adjacency matrix eigenvectors will be processed. Otherwise a discrete laplacian will be used
    :return: biggest element indexes in colums, where each column corresponds to the eigenvector and each column contains nbiggest elements
    """
    if Type == 'Adj':
        absolute = np.absolute(MG.adj_eigenvects)
    else:
        absolute = np.absolute(MG.cond_eigenvects)
    eigbiggest = np.argsort(absolute, axis=0)[-nbiggest:, :]

    for i in range(0, eigbiggest.shape[1]):
        for index, value in reversed(zip(eigbiggest[:, i], absolute[eigbiggest[:, i], i])):
            print index, value, MG.get_descriptor_for_index(index)
        print '<==================>'
        print 'Index \t value \t Descriptor'

    return eigbiggest


def columnSort():
    """
    Computes the degree of adjacency for each node and writes it out in a human-readable format to a file specified in
    Dumps.Adj_degree
    """
    summat = MG.Ajacency_Matrix.sum(axis=1)
    sorts = summat.argsort(axis=0)

    outf = file(Dumps.Adj_degree,'w')

    for i in reversed(range(0,sorts.shape[0])):

        Stri = str(MG.MatrixNumber2NodeID[sorts[i,0]]) + '\t' + str(summat[sorts[i,0],0]) + '\t'+str(MG.get_descriptor_for_index(sorts[i,0]))+'\n'
        outf.write(Stri)

    outf.close()


def compute_sample_circulation_intensity_minimal(Sample, epsilon=1e-10):
    """
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation.

    The informativities are calculated only by using the uniprot proteins as the source and extraction points.
    This reduces the number of interations from  25k to 5 and the number of LU decompositions in a similar manner
    
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    :return: the informativity array
    """

    InformativityArray = np.zeros(( MG.Conductance_Matrix.shape[0], 1))
    Solver = cholesky(csc_matrix(MG.Conductance_Matrix), epsilon)

    print 1, MG.time()
    print Sample
    print len(Sample)

    for i,j in combinations(Sample,2):
        J = np.zeros((MG.Conductance_Matrix.shape[0], 1))
        J[i, 0] = 1.0
        J[j,0] = -1.0
        V = Solver(J)
        Current = get_Current_all(V, J)
        InformativityArray += Current
    return InformativityArray


def compute_sample_circulation_intensity(Sample, epsilon=1e-10, array_v=''):
    """
    Attention, this method shoyld never be called. Normally all the major IO operations on information arrays should be performed with the compute and store
    circulation method

    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    :param array_v: Array version, allowing it to vbe differentiated from other arrays on storage
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    """
    InformativityArray = compute_sample_circulation_intensity_minimal(Sample, epsilon)
    name = Dumps.InfoArray+str(array_v)+'.dump'
    pickle.dump(InformativityArray, file(name, 'w'))
    return InformativityArray


def Compute_circulation_intensity(epsilon=1e-10):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation

    :warning: this method performs the computation on the whole set Uniprots, and thus it's convergence time as a single process is overa month.
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    
    '''

    return compute_sample_circulation_intensity(MG.Uniprot_Mat_idxs, epsilon, '_Full')


def Compute_random_sample(sample_size, iterations, epsilon=1e-10, name_version=''):
    '''
    Performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for
    more rapid computation of information circulation
    Unlike them Method above, if the sample size is under 170, it converges in an hour instead of months.

    :warning: it's computations are not reliable, since it computes a full informativity within a restricted set of proteins.
    :param sample_size: Size of elements to randomly pull out of the uniprot list
    :param iterations: Number of random pulls to perform from the uniprot list to get all the interesing thigns
    :param name_version: suffix to add to the name before the storage to distinguish from the previous iterations.
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float

    '''
    re_UP = copy.copy(MG.Uniprot_Mat_idxs)
    for i in range(0, iterations):
        random.shuffle(re_UP)
        re_UP = re_UP[:sample_size]
        compute_sample_circulation_intensity(re_UP, epsilon, name_version+str(i))


def Compute_truly_random_sample(rounds, iterations, epsilon, name_version=''):
    """
    Performs a computation of information circulation for pairs of randomly selected pairs of nodes.
    Unlike the previous method, it doesn't have a group biais; and the whole matrix is supressed for n repetitions of
    information transfer. Only the nodes that transport any information in addition to their own (i.e.) non-isolated
    and non-leaf nodes. The result is returned and dumped as a new array that is name-versionned with the name_version
    argument.

    Historically this method was used to calibrate the information circulation and rapidly find proteins that were
    transporting the most information

    :param rounds: times random pairs will be selected
    :param iterations: full rounds of informativity computation, each time informativity being added. More or less redoundant with rounds: the total number of iterations is iterations*rounds
    :param epsilon: fudge, preventing overflow
    :param name_version: calculation-specific string allowing to differnetiate the dump from the rest
    """
    Solver = cholesky(csc_matrix(MG.Conductance_Matrix), epsilon)

    for i in range(0, iterations):
        InformativityArray = np.zeros((MG.Conductance_Matrix.shape[0], 1))
        print 'iteration', i, MG.time()

        List_of_pairs = []
        for _ in repeat(None,rounds):
            L = MG.Uniprot_Mat_idxs.copy()
            random.shuffle(L)
            List_of_pairs += zip(L[:len(L)/2], L[len(L)/2:])

        j = 0
        for pair in List_of_pairs:
            j += 1
            if j % 100 == 99:
                print '*'
            J = np.zeros((MG.Conductance_Matrix.shape[0], 1))
            J[pair[0], 0], J[pair[1], 0]  = (1.0, -1.0)
            V = Solver(J)
            Current = get_Current_all(V, J)
            InformativityArray += Current

        CorrInf = np.zeros((MG.Conductance_Matrix.shape[0], 1))
        CorrInf[:] = rounds # names how many times each element was called to communicate with a different one.
        InformativityArray = InformativityArray - CorrInf
        name = Dumps.InfoArray+str(name_version)+str(i)+'.dump'
        Fle = file(name,'w')
        pickle.dump(InformativityArray,Fle)
        Fle.close()
        

def Analyze_relations(Current, number, AdditionalInfos=None):
    """
    Just a way to print out the nodes making pass the most of the current.

    Historically a reporting function, not essential for the matrix manipulation. Should be moved elsewhere.

    :param Current: current running through the graph laplacian of a current system.
    :type Current: numpy array
      
    :param number: the number of the most important relations one wants to observe
    :type number: integer

    :param AdditionalInfos: controls if additional informations will be added to the output. Defaults to None

    """
    # TODO: move this function elsewhere; it is a formatting-for output, not essential for processing function
    # TODO: replace buffered string concatenation by ','.join(arglist)

    writer = file('Current_Analysis.csv','w')
    infos = np.absolute(Current)
    initLine = 'ID'+'\t'+'Mean Informativity'+'\t'+'Type'+'\t'+'Standard deviation'+'\t'+'Pessimistic Info estimation'
    initLine += '\t'+'optimistic Info estimation'+'\t'+'Protein Aboundace'+'\t'+'Essential for Overington'+'\t'+'GO Name'
    initLine += '\t'+'GO ID'+'\t'+'Score'+'\t'+'displayName'+'\t'+'localization'+'\n'
    writer.write(initLine)
    if AdditionalInfos is None:
        for i in range(0,min(number,len(infos))):
            locmax = np.argmax(infos)
            Buffer = ''
            Buffer += str(MG.MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
            Buffer += str(MG.ID2Type[MG.MatrixNumber2NodeID[locmax]])+'\t'
            infos[locmax] = 0.0
            writer.write(Buffer)
    else:
        for i in range(0,min(number,len(infos))):
            locmax = np.argmax(infos)
            Buffer = ''
            Buffer += str(MG.MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
            Buffer += str(MG.ID2Type[MG.MatrixNumber2NodeID[locmax]])+'\t'
            for elt in AdditionalInfos[locmax,:].tolist():
                if str(elt) != str(0.0):
                    Buffer += str(elt)+'\t'
                else:
                    Buffer += '\t'
            Buffer += str(MG.ID2displayName[MG.MatrixNumber2NodeID[locmax]])
            if MG.MatrixNumber2NodeID[locmax] in MG.ID2Localization.keys():
                Buffer += '\t'+str(MG.ID2Localization[MG.MatrixNumber2NodeID[locmax]])
            Buffer += '\n'
            infos[locmax] = 0.0
            writer.write(Buffer)
    writer.close()


def Info_circulation_for_Single_Node(source_MatrixID, sinks_MatrixIDs, epsilon=1e-10):
    '''
    Computes the information circulation for a single node, according to a slightly improved method compared to the one described by Missiuro

    :param source_MatrixID: number of the row/line of the source node within the conductance matrix corresponding to the source
    :type source_MatrixID: int
    
    :param sinks_MatrixIDs: list of numbers of the rows/lines within the conductance matrix corresponding to sinks
    :type sinks_MatrixIDs: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    '''
    Solver = cholesky(csc_matrix(MG.Conductance_Matrix), epsilon)
    Cummulative_Informativity = np.zeros((MG.Conductance_Matrix.shape[0], 1))
    for sink in sinks_MatrixIDs:
        J = np.zeros((MG.Conductance_Matrix.shape[0], 1))
        J[source_MatrixID, 0] = 1.0
        J[sink, 0] = -1.0
        Voltage = Solver(J)
        current = get_Current_all(Voltage, J)
        Cummulative_Informativity += current
    return Cummulative_Informativity


def get_Current_all(Voltages, J):
    '''
    Recovers the current for all the nodes
    
    :param J: Conductance matrix
    :type J: numpy array
    
    :param Voltages: Informativity gradient in each node obtained by the solution of the matrix equation ConductanceMatrix*Voltage = J
    :type Voltages: numpy array

    :return: Currents throug each node.
    :rtype: csc_matrix
    '''
    # This one seems to be pretty f#-up, likely because of the non-connexity of a large portion of uniprots.
    diag_Voltages = lil_matrix(diags(Voltages.T.tolist()[0], 0))
    Corr_Conductance_Matrix = MG.Conductance_Matrix - lil_matrix(diags(MG.Conductance_Matrix.diagonal(), 0))
    sm = diag_Voltages.dot(Corr_Conductance_Matrix) - Corr_Conductance_Matrix.dot(diag_Voltages)
    sign = sm.sign()
    absol = sm.multiply(sign)
    Currents = (np.array(absol.sum(axis=0).T) + np.absolute(J)) / 2.0
    return Currents


def check_Silverality(sample_size, iterations):
    """
    Checks how rapidly the information passing through the neighbouring proteins decreases
    with the distance between them

    :param sample_size: number of roots for the Silverality verification
    :param iterations: how many passes for the Silverality verification
    :return:
    """
    from scipy.sparse.csgraph import dijkstra
    cumulative = []
    rei_UP = copy.copy(MG.Uniprot_Mat_idxs)
    for i in range(0, iterations):
        random.shuffle(rei_UP)
        re_UP = rei_UP[:sample_size+1]
        source_MatrixID = re_UP[0]
        sinks_MatrixIDs = re_UP[0:]
        Currents = Info_circulation_for_Single_Node( source_MatrixID, sinks_MatrixIDs, epsilon=1e-10)
        random.shuffle(rei_UP)
        targets = rei_UP[:sample_size/10]
        distances = dijkstra(MG.Ajacency_Matrix, indices=source_MatrixID, unweighted=True)
        for index in targets:
            cumulative.append((distances[index], Currents[index,0]))
        for i in range(0,sample_size/5):
            index = np.argmin(distances)
            print distances
            cumulative.append((distances[index], Currents[index,0]))
            distances[index] = 100
    pickle.dump(cumulative, file(Dumps.Silverality, 'w'))
    return cumulative


def analyze_Silverality():
    """
    Plots the results of the Silverality Analysis, after undumping them

    """
    Silverality_List = pickle.load(file(Dumps.Silverality,'r'))
    SuperDict = {}
    for distance, value in Silverality_List:
        if distance not in SuperDict.keys():
            SuperDict[distance] = []
        SuperDict[distance].append(value)
    
    for key, val in SuperDict.iteritems():
        print key, val
    
    ArrayDict = {}
    StatsDict = {}
    x = np.zeros((10,1))
    y1 = np.zeros((10,1))
    y2 = np.zeros((10,1))
    y3 = np.zeros((10,1))
    i = 0
    Srtd = sorted(SuperDict.iteritems(), key=operator.itemgetter(0))
    for key, val in Srtd:
        if i > 9:
            break
        ArrayDict[key] = np.array(val)
        StatsDict[key] = (np.mean(ArrayDict[key]),np.std(ArrayDict[key]))
        x[i,0] = key
        y1[i,0] = StatsDict[key][0]
        y2[i,0] = StatsDict[key][0]+StatsDict[key][1]
        y3[i,0] = StatsDict[key][0]-StatsDict[key][1]
        i += 1

    pylab.plot(x, y1, '-k', label='mean')
    pylab.plot(x, y2, '-r', label='mean+std')
    pylab.plot(x, y3, '-b', label='mean-std')
    pylab.legend(loc='upper right')
    pylab.show()



def Compute_and_Store_circulation(handle, List_of_UPs, mongo_db_collection):
    """
    Computes information circulation for a set of uniprots, then pickles the result stores them in a mongoDB

    :param handle: name of the entity to store under the ID
    :param List_of_UPs: list of Uniprot Nodes on which we are willing to perform the retrieval and storage
    :param mongo_db_collection:  pymongo mongoDB client attachement
    :return: Informativity array that was stored
    """
    print 'entering computation for: ', handle, 'with', len(List_of_UPs), 'UPs'
    re_UP_List = []
    for elt in List_of_UPs:
        if elt in MG.NodeID2MatrixNumber.keys():
            re_UP_List.append(MG.NodeID2MatrixNumber[elt])

    Informativity_Array = compute_sample_circulation_intensity_minimal(re_UP_List)
    UPs_Set = pickle.dumps(set(List_of_UPs))
    post = {'GO_ID': handle, 'UP_Set':UPs_Set, 'Info_Array':pickle.dumps(Informativity_Array)}
    mongo_db_collection.insert(post)
    return Informativity_Array


def Compute_ponderated_info_circulation(UPs_2Binding_Affs):
    """
    Should compute the information circulation wehre the flowing information is ponderated by the importance a given
    protein has for our analysis.

    :param UPs_2Binding_Affs: Dictionary linking UniprotNode_IDs to their relative importances for the ponderated information circ. calculation
    :raise NotImplementedError: It is not implemented yet
    """
    raise NotImplementedError


if __name__ == "__main__":
    # characterise_eigenvects(100, Type = 'Cond')
    # checkMatrix()
    processEigenVectors(15,Type = 'Cond')
    # columnSort()
    # lst = [22811, 18147, 13023]
    # print compute_sample_circulation_intensity_minimal(lst)
    # check_Silverality(10, 100)
    # analyze_Silverality()
    pass


# TODO: for a double node and a single node computation, add internode voltage in addition to the information flow.
    # the formula is actually J_convolve_Voltages_times_J.T

# TODO: Compute the voltages during the computation to determine what would be the current between two nodes at a given potential difference.
# TODO: Compute the currents at a given voltage for each GO term pair.
# TODO: Create GO_Term - specific proximity matrix, and store it with the others for the later clustering

# TODO: compute the purities of main action: 3 most important contibutions to GO terms, then pull them into a correct position and comprare
    # - Importance
    # - purity
    # - possible clustering => use the BEA algo from the Part 1 of the internship

# TODO: compute the whole ensemble of the nodes perturbed by the protein
        # * segment them with the voltage-based interaction
        # * retrieve segments
        # * for each segment compute the most critical protein (hidden variable)
        # * compute the annotation ponderated by the importance for all of the proteins significantly affectd by the conductance calculation
        # * (don't forget to divide by n each node), where n is the
        # * compare the significance of the protein affected to the information circulation in the whole interactome

# TODO: implementation while using Ehit interactions only
