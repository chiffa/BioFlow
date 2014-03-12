'''
Created on Jul 11, 2013
@author: andrei

This module contains all the routines that are respojnsible for pulling 
the matrixes out of the neo4j graph and processing them

'''
from PolyPharma.configs import Dumps
import copy
from scipy.sparse import lil_matrix
# from scipy.sparse.linalg import eigsh
# import itertools
# from time import time
import pickle
import numpy as np
import operator
# from os import listdir
import random
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky
from scipy.sparse import csc_matrix
from scipy.sparse import diags
ID2Aboundances = {}
import pylab
# import math
from Matrix_Interactome_DB_interface import MatrixGetter


MG = MatrixGetter(True, False)
MG.fast_load()


# TODO: we need to refactor the filtering system

def _characterize(Objkt):
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
    sbigest_args = np.argsort(norm)[-N_of_eigenvectors_to_characterise:]

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
    '''
    Computes the degree of adjacency for each node and writes it out in a human-readable format to a file specified in
    Dumps.Adj_degree
    '''

    summat = MG.Ajacency_Matrix.sum(axis=1)
    sorts = summat.argsort(axis=0)

    outf = file(Dumps.Adj_degree,'w')

    for i in reversed(range(0,sorts.shape[0])):

        Stri = str(MG.MatrixNumber2NodeID[sorts[i,0]]) + '\t' + str(summat[sorts[i,0],0]) + '\t'+str(MG.get_descriptor_for_index(sorts[i,0]))+'\n'
        outf.write(Stri)

    outf.close()


# def get_voltages(numpy_array, InformativityDict):
#     """
#     Deprecated.
#     This method seems to compute the voltages from several computation
#
#     :param numpy_array: array of voltages on each node resulting from the computation
#     :param InformativityDict:
#     """
#     for i in range(0, len(numpy_array)):
#         InformativityDict[MG.MatrixNumber2NodeID[i]] += numpy_array[i, 0]
#
#
#
# def create_InfoDict():
#     """
#     Deprecated
#     Creates a dict mapping all the DB NodeIDs to zero. Might be a home-brewed alternative to collections.defaultdict
#
#     :return: dict
#     """
#     new_dict={}
#     for val in MG.MatrixNumber2NodeID.itervalues():
#         new_dict[val]=0.0
#     return new_dict


def compute_sample_circulation_intensity_minimal(Sample, epsilon=1e-10):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation.

    The informativities are calculated only by using the uniprot proteins as the source and extraction points.
    This reduces the number of interations from  25k to 5 and the number of LU decompositions in a similar manner
    
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    :return: the informativity array
    
    '''

    InformativityArray = np.zeros(( MG.Conductance_Matrix.shape[0], 1))                                # Database ID to Informativity
    Solver = cholesky(csc_matrix(MG.Conductance_Matrix), epsilon)

    print 1, MG.time()
    print Sample
    print len(Sample)
    for i in range(0,len(Sample)):
        print '\n', 2, i, MG.time()
        for j in range(i,len(Sample)):
            if j%25 == 24:
                print '*',
            J = np.zeros((MG.Conductance_Matrix.shape[0], 1))   #TODO: refactor to use itertools.combinations
            J[Sample[i],0] = 1.0
            J[Sample[j],0] = -1.0
            V = Solver(J)
            Current = get_Current_all(V, J)
            InformativityArray += Current
    return InformativityArray


def compute_sample_circulation_intensity(Sample, epsilon=1e-10, array_v=''):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    :param array_v:
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    
    '''
    # NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    InformativityArray = compute_sample_circulation_intensity_minimal(Sample, epsilon)
    name = 'InfoArray_new'+str(array_v)+'.dump'   # TODO: change the undump behavior here to a MongoDB call
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
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for
    more rapid computation of information circulation

    Unlike them ethod above, if the sample size is under 170, it converges in an hour instead of months.

    :warning: it's computations are not reliable, since it computes a full informativity within a restricted set of proteins.

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
    Performs a computation of information circulation for pairs of randomly selected pairs of nodes

    :param rounds:
    :param iterations:
    :param epsilon:
    :param name_version:
    """
    Solver = cholesky(csc_matrix(MG.Conductance_Matrix),epsilon)

    for i in range(0, iterations):
        InformativityArray = np.zeros((MG.Conductance_Matrix.shape[0], 1))
        print 'iteration', i, MG.time()

        L1, L2 = (MG.Uniprot_Mat_idxs.copy(),MG.Uniprot_Mat_idxs.copy())
        random.shuffle(L1)
        random.shuffle(L2)
        List_of_pairs = zip(L1, L2)

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
        name = 'InfoArray_new'+str(name_version)+str(i)+'.dump'   # TODO: change the undump behavior here to a MongoDB call
        Fle = file(name,'w')
        pickle.dump(InformativityArray,Fle)
        Fle.close()
        

def Analyze_relations(Current, number, AdditionalInfos=None):
    '''
    Just a way to print out the nodes making pass the most of the current
    
    :param Current: current running through the graph laplacian of a current system.
    :type Current: numpy array
      
    :param number: the number of the most important relations one wants to observe
    :type number: integer
    
    '''
    writer = file('Current_Analysis.csv','w')
    infos = np.absolute(Current)
    initLine = 'ID'+'\t'+'Mean Informativity'+'\t'+'Type'+'\t'+'Standard deviation'+'\t'+'Pessimistic Info estimation'
    initLine += '\t'+'optimistic Info estimation'+'\t'+'Protein Aboundace'+'\t'+'Essential for Overington'+'\t'+'GO Name'
    initLine += '\t'+'GO ID'+'\t'+'Score'+'\t'+'displayName'+'\t'+'localization'+'\n'

    writer.write(initLine)
    if AdditionalInfos == None:
        for i in range(0,min(number,len(infos))):
            locmax = np.argmax(infos)
            Buffer = ''
            Buffer += str(MG.MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
            print MG.MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer += str(MG.ID2Type[MG.MatrixNumber2NodeID[locmax]])+'\t'
            print MG.ID2Type[MG.MatrixNumber2NodeID[locmax]], '\t',
            Buffer += str(MG.ID2displayName[MG.MatrixNumber2NodeID[locmax]])+'\n'
            print MG.ID2displayName[MG.MatrixNumber2NodeID[locmax]]
            infos[locmax] = 0.0
            writer.write(Buffer)
    else:
        for i in range(0,min(number,len(infos))):
            locmax = np.argmax(infos)
            Buffer = ''
            Buffer += str(MG.MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
            print MG.MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer += str(MG.ID2Type[MG.MatrixNumber2NodeID[locmax]])+'\t'
            print MG.ID2Type[MG.MatrixNumber2NodeID[locmax]], '\t',
            for elt in AdditionalInfos[locmax,:].tolist():
                if str(elt) != str(0.0):
                    Buffer += str(elt)+'\t'
                    print elt, '\t',
                else:
                    Buffer += '\t'
            Buffer += str(MG.ID2displayName[MG.MatrixNumber2NodeID[locmax]])
            if MG.MatrixNumber2NodeID[locmax] in MG.ID2Localization.keys():
                Buffer += '\t'+str(MG.ID2Localization[MG.MatrixNumber2NodeID[locmax]])
            Buffer += '\n'
            print MG.ID2displayName[MG.MatrixNumber2NodeID[locmax]]
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
    # This one seems to be pretty fucked-up, likely because of the non-connexity of a large portion of uniprots.
    diag_Voltages = lil_matrix(diags(Voltages.T.tolist()[0], 0))
    Corr_Conductance_Matrix = MG.Conductance_Matrix - lil_matrix(diags(MG.Conductance_Matrix.diagonal(), 0))
    sm = diag_Voltages.dot(Corr_Conductance_Matrix) - Corr_Conductance_Matrix.dot(diag_Voltages)
    sign = sm.sign()
    absol = sm.multiply(sign)
    Currents = (np.array(absol.sum(axis=0).T) + np.absolute(J)) / 2.0
    return Currents

    
# def stats_over_random_info_circ_samples(UniProtAttachement=True):
#     from PolyPharma.Utils.Prot_Aboundances import ID2Aboundances
#     UPNode_IDs_2Proteins_IDs_List = MG.Uniprot_attachments
#     NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))   # TODO: change the undump behavior here
#     DicList=[]
#     FnameList=[]
#     Filenames = listdir('.')
#     for filename in Filenames:    # TODO: change the undump behavior here
#         if 'InfoArray_new' in filename:
#             FnameList.append(filename)
#     for fname in FnameList:
#         DicList.append(pickle.load(file(fname,'r')))
#     Stats_Mat=np.concatenate(tuple(DicList),axis=1)
#     MeanInfos=np.mean(Stats_Mat,axis=1).reshape((Stats_Mat.shape[0],1))
#     STDInfos=np.std(Stats_Mat,axis=1).reshape((Stats_Mat.shape[0],1))
#     # Check if STD>MEANInfos
#     for LineID in range(0,Stats_Mat.shape[0]):
#         if MeanInfos[LineID,0] < STDInfos[LineID,0]:
#             print LineID, MatrixNumber2NodeID[LineID],'|',MeanInfos[LineID,0], STDInfos[LineID,0]
#             print Stats_Mat[LineID,:]
#             print '<==============================>'
#     if UniProtAttachement:
#         # @attention: this is a temporary patch.
#         # TODO: A full implementation would require several iterations of eHiT propagation and then
#         # several more iterations of Reactome.org propagation.
#         for UP_ID in UPNode_IDs_2Proteins_IDs_List.keys():
#             if len(UPNode_IDs_2Proteins_IDs_List[UP_ID])>1:
#                 print UP_ID, 'problem!!!!'
#             else:
#                 Prot_ID=UPNode_IDs_2Proteins_IDs_List[UP_ID][0]
#                 if Prot_ID in NodeID2MatrixNumber.keys():
#                     MeanInfos[NodeID2MatrixNumber[UP_ID],0]+=MeanInfos[NodeID2MatrixNumber[Prot_ID],0]
#                     ID2Localization[UP_ID]=ID2Localization[Prot_ID]
#                     STDInfos[NodeID2MatrixNumber[UP_ID],0]=math.sqrt(STDInfos[NodeID2MatrixNumber[UP_ID],0]**2+STDInfos[NodeID2MatrixNumber[Prot_ID],0]**2)
#     pessimist=MeanInfos-1.97*STDInfos
#     optimist=MeanInfos+1.97*STDInfos
#     aboundances=np.zeros((Stats_Mat.shape[0],1))
#     errcount1=0
#     for key,val in ID2Aboundances.iteritems():
#         if key in NodeID2MatrixNumber.keys():
#             aboundances[NodeID2MatrixNumber[key],0]=val
#         else:
#             errcount1+=1
#     Overingtonicitiy=np.zeros((Stats_Mat.shape[0],1))
#     Overington_IDList=pickle.load(file('IDList.dump','r'))   # TODO: change the undump behavior here
#     errcount2=0
#     finmatrix=pickle.load(file('finmatrix.dump','r'))  # TODO: change the undump behavior here
#     for ID in Overington_IDList:
#         if ID in NodeID2MatrixNumber.keys():
#             Overingtonicitiy[NodeID2MatrixNumber[ID],0]=1.0
#         else:
#             errcount2+=1
#     print 'errcount from stats', errcount1, errcount2
#     print STDInfos.shape, pessimist.shape, optimist.shape, aboundances.shape, Overingtonicitiy.shape, finmatrix.shape
#     Additional=np.concatenate((STDInfos,pessimist,optimist,aboundances,Overingtonicitiy,finmatrix),axis=1)
#     Analyze_relations(MeanInfos, 50000, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Additional)


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

#     array1=np.zeros((len(Silverality_List),1))
#     array2=np.zeros((len(Silverality_List),1))
#     
#     for i in range(0, len(Silverality_List)):
#         array1[i,0]=Silverality_List[i][0]
#         array2[i,0]=Silverality_List[i][1]
#     pylab.plot(array1, array2)
#     pylab.show()
    # TODO: perform statistical analysis


def Compute_and_Store_circulation(handle, List_of_UPs, mongo_db_collection):
    """

    :param handle:
    :param List_of_UPs:
    :param mongo_db_collection:
    :return:
    """
    print 'entering computation for: ', handle, 'with', len(List_of_UPs), 'UPs'
    re_UP_List = []
    for elt in List_of_UPs:
        if elt in MG.NodeID2MatrixNumber.keys():
            re_UP_List.append(MG.NodeID2MatrixNumber[elt])
    Informativity_Array = compute_sample_circulation_intensity_minimal(re_UP_List)
    UPs_Set = pickle.dumps(set(List_of_UPs))
    post = {'GO_ID':handle, 'UP_Set':UPs_Set, 'Info_Array':pickle.dumps(Informativity_Array)}
    mongo_db_collection.insert(post)
    return Informativity_Array


def Compute_ponderated_info_circulation(UPs_2Binding_Affs):
    raise NotImplementedError


if __name__ == "__main__":
    # characterise_eigenvects(100, Type = 'Cond')
    # checkMatrix()
    # processEigenVectors(15,Type = 'Cond')
    # columnSort()
    lst = [22811, 18147, 13023]
    compute_sample_circulation_intensity_minimal(lst)
    check_Silverality(10, 100)
    analyze_Silverality()




# TODO: compute the whole ensemble of the nodes perturbed by the protein
        # segment them with the voltage-based interaction
        # retrieve segments
        # for each segment compute the most critical protein (hidden variable)
        # compute the annotation ponderated by hte importance for all of the proteins significantly affectd by the conductance calculation 
        # (don't forget to divide by n each node), where n is the 
        # compare the significance of the protein affected to the information circulation in the whole interactome




# TODO: compute the ponderated Information circulation: use protein aboundances. The other one is assumed to be inhibition power biasis
        # then extract summary GO terms based on the Uniprot affection data

# TODO: along with Overingtonicity integrate the list of essential genes in human diseases from the PLoS 2011 publication

# DONE: reverse GO_Access: provided the Uniprots find the proteins carrying over the most information
# DONE: mount a PyMongo data store in order to be able to save and retrieve the programming objects easily
#         How is it done: - picket to string
#        Store an object in a collection defined by it's Id and computation number
#        If requested, retrieve by ID or else
#         Index on the GO ID and belonging UNIPROTs (If same set of uniprots, it is the same) => store as sets
#         Pickles of sets with the same elements are always the same

# Importance of complementation of the information with the Reactome.org data with the EHiT data: otherwise the information circulation completely sucks
# Reactome.org: the interactions due to kinases aren't explicitly shown. Instead a broadcasting through the secondary features that perform the modification
# Is needed. Which is completely stupid, because it doesn't show the specific action on the proteins due to the conformation modification. Thus Reactome.org
# is more of a ressource for human experts then for truly machine-learning tasks.

# TODO:
# Compute the voltages during the computation to determine what would be the current between two nodes at a given current. 
# Compute the currents at a given voltage for each GO term pair.
# Create GO_Term - specific proximity matrix, and store it with the others for the later clustering

# TODO: compute the purities of main action: 3 most important contibutions to GO terms, then pull them into a correct position and comprare
# - Importance
# - purity
# - possible clustering => use the BEA algo from the Part 1 of the internship

# Perform_Loading_Routines()
# 
# stats_over_random_info_circ_samples(True)
# 
# check_Silverality(100,100)
# 
# analyze_Silverality()
# Compute_truly_random_sample(2,3,1e-10, '1_')
# 
# stats_over_random_info_circ_samples(True)
# 
# this is going to last for a while. Now we need to get it split among several nodes
#
# <================================================>
#
# checkMatrix()
# 
# get_eigenvect_Stats()
# 
# processEigenVectors()
# 
# mass_Calibrate(6,10,True)
#  
# treat_Calibration()
# 
# columnSort()


# DONE: remake the sampling so it is efficiently 170**2/2 one to one randomly chosen pairs that are calculated, and not the whole 170 ensemble, so that the 
# Informativities actually follow a gaussian distribution


# TODO: There might be an error in the module responsible for linkage between the uniprots and the accession numbers: for instance the 20253 has an annotation with an Acnum, but
# has no Uniprot attached to it within the database

# TODO: create GO and Pathway Structure access
# Calibrate the values so that after ~ 3 transitions the correlation vanishes on average (Follow Pamela Silver Approach) => this is actually the cumulated perturbation of
# two targets that shoudl vanish totally 

# TODO: implementation while using Ehit interactions only


