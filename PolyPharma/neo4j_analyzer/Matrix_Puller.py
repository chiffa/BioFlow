'''
Created on Jul 11, 2013
@author: andrei

This module contains all the routines that are respojnsible for pulling 
the matrixes out of the neo4j graph and processing them

'''
from PolyPharma.configs import IDFilter
import copy
from scipy.sparse import lil_matrix
# from scipy.sparse.linalg import eigsh
# import itertools
from time import time
import pickle
import numpy as np
import operator
from os import listdir
import random
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky
from scipy.sparse import csc_matrix
from scipy.sparse import diags
ID2Aboundances = {}
import pylab
import math
from Matrix_Interactome_DB_interface import MatrixGetter

MG = MatrixGetter(True,False)
MG.fast_load()


# TODO: we need to refactor the filtering system
def characterise_eigenvects():
    """
    Recovers statistics over eigenvectors in order to remove excessively connected nodes.

    :return: Index counting how many times each element has occured in the 100 biggest eigenvectors, then sorts them,
             Index characterising each element from the matrix
    """

    # SuperIndex = {}
    # CounterIndex = {}
    # eigenVectsList = np.split(MG.adj_eigenvects,MG.adj_eigenvects.shape[1], axis = 1)

    # TODO: We can do it waay more easily with numpy
    norm = np.linalg.norm(MG.adj_eigenvects, axis=1)
    sbigest_args = np.argsort(norm)[-100:]

    for arg in sbigest_args:
        print arg, norm[arg], MG.get_descriptor_for_index(arg)

    print sbigest_args

    # for eigenvect in eigenVectsList:
    #     normalized = np.multiply(eigenvect, eigenvect)
    #     for k in range(0,10):
    #         index = np.argmax(normalized, 0)
    #         if MG.MatrixNumber2NodeID[index] not in CounterIndex.keys():
    #             CounterIndex[MG.MatrixNumber2NodeID[index]] = 0
    #             SuperIndex[MG.MatrixNumber2NodeID[index]] = MG.get_descriptor_for_index(index)
    #         CounterIndex[MG.MatrixNumber2NodeID[index]] += normalized[index]
    #         normalized[index] = 0
    # srtd = sorted(CounterIndex.iteritems(), key = operator.itemgetter(1), reverse = True)
    # print IDFilter
    # for key, val in srtd[:100]:
    #     print val, key, SuperIndex[key], key in IDFilter
    #     if key in IDFilter:
    #         print 'error on key: ', key

    # return CounterIndex, SuperIndex


def UniprotCalibrate(rounds,depth, filename, Rdom):
    '''
    Checks if the decrease corresponds on average to the value predicted for natural networks by 
    P. Silver in E.Coli. One single propagation iteration
    # NOTICE: it might be better to perform the information collection only for the other uniprot- proteins
    '''
    init = time()
    pickleDump2 = file('pickleDump2.dump','r')  # TODO: change the undump behavior here
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    pickleDump3 = file('pickleDump3.dump','r')  # TODO: change the undump behavior here
    ValueMatrix = pickle.load(pickleDump3)
    print 'unpickled in:', time()-init
    Finale = MainUprotCalLoop(Uniprots, ValueMatrix, NodeID2MatrixNumber, rounds, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, depth,Rdom)
    write = file(filename, 'w')
    pickle.dump(Finale, write)
    write.close()
    print 'round completed in:', time()-init      


def MainUprotCalLoop(Uniprots, ValueMatrix, NodeID2MatrixNumber, rounds, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, depth, Rdom):
    iterations=1
    ReUniprots=copy.copy(Uniprots)
    Finale=[]
    if Rdom:
        iterations=3
    for i in range(0,iterations):
        Subfinale=[]
        if Rdom:
            Uniprots=random.sample(Uniprots, 200)
        for ID in Uniprots:
            Vector=np.zeros((ValueMatrix.shape[1],1))
            Vector[NodeID2MatrixNumber[ID]]=1.0
            ForbidList=[]
            for i in range(0,rounds):
                for k in Vector.nonzero():
                    ForbidList.append(k)
                Vector=ValueMatrix*Vector
                for k in ForbidList:
                    Vector[k]=0.0
            LocalSample={ID:MG.get_descriptor_for_index(NodeID2MatrixNumber[ID])}
            nzIndexes = Vector.nonzero()[0]
            redepth=min(depth,len(nzIndexes))
            suplist=random.sample(nzIndexes,redepth)
            for index in suplist:
                index=int(index)
                value=0
                if np.linalg.norm(Vector,1) < 10e-3:
                    print 'error', ID, np.linalg.norm(Vector,1)
                    value='infinity'
                else:
                    value=Vector[index]/np.linalg.norm(Vector,1)
                LocalSample[MatrixNumber2NodeID[index]]=(Vector[index], value, MG.get_descriptor_for_index(index))
            Subfinale.append(LocalSample)
        Finale.append(Subfinale)
    return Finale


def mass_Calibrate(maxrange,depth,Rdom=False):
    '''
    Performs several rounds of calibration, actually recomputing the figure presented by P. Silver
    '''
    for i in range(2, maxrange+1):
        filename=str('calibrate'+str(i)+'.dump')
        UniprotCalibrate(i,depth,filename,Rdom)


def treat_Calibration():
    '''
    Performs the analysis of calibration
    '''
    DicList={}
    FnameList=[]
    Filenames = listdir('.')
    for filename in Filenames:
        if 'calibrate' in filename:
            FnameList.append(filename)
    for fname in FnameList:               # TODO: change the undump behavior here
        DicList[fname.split('.')[0][-1]]=(pickle.load(file(fname,'r')))
    memory={}    
    for key in DicList.keys():
        memory[key]=[]
        print key, 'iterated!'
        for sublist in DicList[key]:
            average1=0
            average2=0
            count=0 
            for Dict in sublist:
                print Dict
                for subkey, subval in Dict.iteritems():
                    if subval[0]!='UNIPROT':
                        if int(subval[0])<0:
                            print Dict
                        count+=1
                        average1+=int(subval[0])
                        average2+=int(subval[1])
            average1=float(average1)/float(count)
            average2=float(average2)/float(count)
            memory[key].append((average1,average2))
    srtd=sorted(memory.iteritems(), key=operator.itemgetter(0))
    for key, val in srtd:
        print key, val[0][0], val[1][0], val[2][0], '|', val[0][1], val[1][1], val[2][1]
            

def checkMatrix():
    pickleDump3=file('pickleDump3.dump','r')  # TODO: change the undump behavior here
    ValueMatrix=pickle.load(pickleDump3)
    NzeroList=ValueMatrix.nonzero()
    faultyList=[]
    for i in range(0,len(NzeroList[0])):
        index=(int(NzeroList[0][i]),int(NzeroList[1][i]))
        value=ValueMatrix[index[0],index[1]]
        if value < 0.0 or value > 1.0:
            faultyList.append(index)
    print len(faultyList)
    print len(NzeroList[0])
#     rsample=random.sample(range(0,len(NzeroList[0])), 10)
#     for i in rsample:
#         index=(int(NzeroList[0][i]),int(NzeroList[1][i]))
#         value=ValueMatrix[index[0],index[1]]
#         print index, value


def processEigenVectors():
    init=time()
    pickleDump=file('pickleDump.dump','r')  # TODO: change the undump behavior here
    eigenvals, eigenvects=pickle.load(pickleDump)
    pickleDump2=file('pickleDump2.dump','r')  # TODO: change the undump behavior here
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    eigenVectsList=[]
    SuperIndex=[]
    print 'depicked',time()-init
    for i in range(0,eigenvects.shape[1]):
        eigenVectsList.append(eigenvects[:,i])
    i=0
    for eigenvect in eigenVectsList:
        i+=1
        normalized=np.multiply(eigenvect,eigenvect)
        indexes={}
        for k in range(0,10):
            index=np.argmax(normalized, 0)
            indexes[MatrixNumber2NodeID[index]]=(normalized[index], MG.get_descriptor_for_index(index))
            normalized[index]=0
        SuperIndex.append(indexes)
    for subIndex in SuperIndex:
        for key in subIndex.keys():
            print '\t', key, subIndex[key]
        print '<=========================>'
    return SuperIndex


def columnSort():
    '''
    np.sum is broken for sparse matrixes. does the same thing with option axis=0
    '''
    pickleDump3=file('pickleDump3.dump','r')  # TODO: change the undump behavior here
    ValueMatrix=pickle.load(pickleDump3)
    SupportDict={}
    IndexDict={}
    nz=ValueMatrix.nonzero()
    pickleDump2=file('pickleDump2.dump','r')  # TODO: change the undump behavior here
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(pickleDump2)
    print 'ended imports'
    for i in range(0, len(nz[0])):
        if nz[0][i] not in SupportDict.keys():
            SupportDict[nz[0][i]]=[]
        SupportDict[nz[0][i]].append(nz[1][i])
    print 'ended loading indexes'
    for key in SupportDict.keys():
        if len(SupportDict[key])>1:
            IndexDict[key]=0
            for val in SupportDict[key]:
                IndexDict[key]+=float(ValueMatrix[key,val])
                
    srtd = sorted( IndexDict.iteritems(), key=operator.itemgetter(1),reverse=True )
    
    outf = file('columnSum.csv','w')
    
    for elt in srtd:
        Stri=str(str(MatrixNumber2NodeID[elt[0]])+'\t'+str(elt[1])+'\n')
        print  MatrixNumber2NodeID[elt[0]], elt[1]
        outf.write(Stri)
    outf.close()


def get_voltages(numpy_array, MatrixNumber2NodeID, InformativityDict):
    for i in range(0, len(numpy_array)):
        InformativityDict[MatrixNumber2NodeID[i]]+=numpy_array[i,0]
    return


def create_InfoDict(MatrixNumber2NodeID):
    new_dict={}
    for val in MatrixNumber2NodeID.values():
        new_dict[val]=0.0
    return new_dict


def compute_sample_circulation_intensity_minimal(Sample, epsilon=1e-10):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    
    '''
    init=time()
    conductance_Matrix=pickle.load(file('pickleDump4.dump','r'))   # TODO: change the undump behavior here
    # NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    InformativityArray=np.zeros((conductance_Matrix.shape[0],1))                                # Database ID to Informativity
    Solver=cholesky(csc_matrix(conductance_Matrix),epsilon)
    # The informativities are calculated only by using the uniprot proteins as the source and extraction points.
    # This reduces the number of interations from  25k to 5 and the number of LU decompositions in a similar manner
    print 1, time()-init
    init=time()
    print Sample
    print len(Sample)
    for i in range(0,len(Sample)):
        print '\n', 2, i, time()-init
        init=time()
        for j in range(i,len(Sample)):
            if j%25==24:
                print '*',
            J=np.zeros((conductance_Matrix.shape[0],1))
            J[Sample[i],0]=1.0
            J[Sample[j],0]=-1.0
            V=Solver(J)
            Current=get_Current_all(conductance_Matrix,V,J)
            InformativityArray+=Current
    return InformativityArray


def compute_sample_circulation_intensity(Sample, epsilon=1e-10, array_v=''):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    :param Sample: List of row/column numbers that are to be used in the computation (~170 for a computation converging in less then an hour)
    :type Sample: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    
    '''
    # NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))
    InformativityArray=compute_sample_circulation_intensity_minimal(Sample, epsilon)
    name='InfoArray_new'+str(array_v)+'.dump'   # TODO: change the undump behavior here
    pickle.dump(InformativityArray,file(name,'w'))
    return InformativityArray


def Compute_circulation_intensity(epsilon=1e-10):
    '''
    performs the information circulation calculation in agreement with the publication by Misiuro et al, but with algorithm slightly modified for 
    more rapid computation of information circulation
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    
    '''
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))   # TODO: change the undump behavior here
    re_UP=[]
    for SP in Uniprots:
        re_UP.append(NodeID2MatrixNumber[SP])
    return compute_sample_circulation_intensity(re_UP,epsilon,'_Full')


def Compute_random_sample(sample_size, iterations, epsilon=1e-10, name_version=''):
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))   # TODO: change the undump behavior here
    re_UP=[]
    for SP in Uniprots:
        re_UP.append(NodeID2MatrixNumber[SP])
    for i in range(0,iterations):
        random.shuffle(re_UP)
        re_UP=re_UP[:sample_size]
        compute_sample_circulation_intensity(re_UP,epsilon,name_version+str(i))


def Compute_truly_random_sample(rounds,iterations,epsilon,name_version=''):
    conductance_Matrix=pickle.load(file('pickleDump4.dump','r'))   # TODO: change the undump behavior here
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))   # TODO: change the undump behavior here
    Solver=cholesky(csc_matrix(conductance_Matrix),epsilon)
    re_UP=[]
    for SP in Uniprots:
        re_UP.append(NodeID2MatrixNumber[SP])
    init=time()
    for i in range(0,iterations):
        InformativityArray=np.zeros((conductance_Matrix.shape[0],1))
        print 'iteration', i, time()-init
        init=time()
        List_of_pairs=[]
        for j in range(0,rounds):
            Lst=copy.copy(re_UP)
            random.shuffle(Lst)
            while len(Lst)>2:
                elt1=Lst.pop()
                elt2=Lst.pop()
                List_of_pairs.append((elt1,elt2))
        j=0
        init2=time()
        for pair in List_of_pairs:
            j+=1
            if j%100==99:
                print '*', time()-init2,
                init2=time()
            J=np.zeros((conductance_Matrix.shape[0],1))
            J[pair[0],0]=1.0
            J[pair[1],0]=-1.0
            V=Solver(J)
            Current=get_Current_all(conductance_Matrix,V,J)
            InformativityArray+=Current
        CorrInf=np.zeros((conductance_Matrix.shape[0],1))
        CorrInf[:]=rounds
        InformativityArray=InformativityArray-CorrInf
        name='InfoArray_new'+str(name_version)+str(i)+'.dump'   # TODO: change the undump behavior here
        Fle=file(name,'w')
        pickle.dump(InformativityArray,Fle)
        Fle.close()
        

def Analyze_relations(Current, number,MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization,AdditionalInfos=None):
    '''
    Just a way to print out the nodes making pass the most of the current
    
    :param Current: the current
    :type Current: numpy array
      
    :param number: the number of the most important relations one wants to observe
    :type number: integer
    
    '''
    writer=file('Current_Analysis.csv','w')
    infos=np.absolute(Current)
    initLine='ID'+'\t'+'Mean Informativity'+'\t'+'Type'+'\t'+'Standard deviation'+'\t'+'Pessimistic Info estimation'+'\t'+'optimistic Info estimation'+'\t'+'Protein Aboundace'+'\t'+'Essential for Overington'+'\t'+'GO Name'+'\t'+'GO ID'+'\t'+'Score'+'\t'+'displayName'+'\t'+'localization'+'\n'
    writer.write(initLine)
    if AdditionalInfos==None:
        for i in range(0,min(number,len(infos))):
            locmax=np.argmax(infos)
            Buffer=''
            Buffer+=str(MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
#             print MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer+=str(ID2Type[MatrixNumber2NodeID[locmax]])+'\t'
#             print ID2Type[MatrixNumber2NodeID[locmax]], '\t',
            Buffer+=str(ID2displayName[MatrixNumber2NodeID[locmax]])+'\n'
#             print ID2displayName[MatrixNumber2NodeID[locmax]]
            infos[locmax]=0.0
            writer.write(Buffer)
    else:
        for i in range(0,min(number,len(infos))):
            locmax=np.argmax(infos)
            Buffer=''
            Buffer+=str(MatrixNumber2NodeID[locmax])+'\t'+str(float(infos[locmax]))+'\t'
#             print MatrixNumber2NodeID[locmax], float(infos[locmax]),'\t',
            Buffer+=str(ID2Type[MatrixNumber2NodeID[locmax]])+'\t'
#             print ID2Type[MatrixNumber2NodeID[locmax]], '\t',
            for elt in AdditionalInfos[locmax,:].tolist():
                if str(elt)!=str(0.0):
                    Buffer+=str(elt)+'\t'
    #                 print elt, '\t',
                else:
                    Buffer+='\t'
            Buffer+=str(ID2displayName[MatrixNumber2NodeID[locmax]])
            if MatrixNumber2NodeID[locmax] in ID2Localization.keys():
                Buffer+='\t'+str(ID2Localization[MatrixNumber2NodeID[locmax]])
            Buffer+='\n'
#             print ID2displayName[MatrixNumber2NodeID[locmax]]
            infos[locmax]=0.0
            writer.write(Buffer)
    writer.close()


def Info_circulation_for_Single_Node(conductance_Matrix,source_MatrixID,sinks_MatrixIDs,epsilon=1e-10):
    '''
    Computes the information circulation for a single node, according to a slightly improved method compared to the one described by Missiuro

    :param conductance_Matrix: Conductance matrix
    :type conductance_Matrix: scipy.sparse lil_matrix
    
    :param source_MatrixID: number of the row/line of the source node within the conductance matrix corresponding to the source
    :type source_MatrixID: int
    
    :param sinks_MatrixIDs: list of numbers of the rows/lines within the conductance matrix corresponding to sinks
    :type sinks_MatrixIDs: list of ints
    
    :param epsilon: corrective factor for the calculation of Cholesky matrix. defaults to 1e-10
    :type epsilon: float
    '''
    Solver=cholesky(csc_matrix(conductance_Matrix),epsilon)
    Cummulative_Informativity=np.zeros((conductance_Matrix.shape[0],1))
    for sink in sinks_MatrixIDs:
        J=np.zeros((conductance_Matrix.shape[0],1))
        J[source_MatrixID,0]=1.0
        J[sink,0]=-1.0
        Voltage=Solver(J)
        current=get_Current_all(conductance_Matrix,Voltage,J)
        Cummulative_Informativity+=current
    return Cummulative_Informativity


def get_Current_all(Conductance_Matrix, Voltages, J):
    '''
    Recovers the current for all the nodes
    
    @param Conductance_Matrix: Conductance matrix
    @type Conducntace_Matrix: scipy.sparse lil_matrix
    
    @param Voltages: Informativity gradient in each node obtained by the solution of the matrix equation ConductanceMatrix*Voltage = J
    @type Voltages: numpy array
    
    '''
    diag_Voltages=lil_matrix(diags(Voltages.T.tolist()[0],0))
    Corr_Conductance_Matrix=Conductance_Matrix-lil_matrix(diags(Conductance_Matrix.diagonal(),0))
    Currents=(np.absolute(diag_Voltages*Corr_Conductance_Matrix-Corr_Conductance_Matrix*diag_Voltages).sum(axis=0).T+np.absolute(J))/2.0
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


def check_Silverality(sample_size,iterations):
    '''
    Checks how rapidly the information passing through the neighbouring proteins decreases
    with the distance between them
    '''
    from scipy.sparse.csgraph import dijkstra
    conductance_Matrix=pickle.load(file('pickleDump4.dump','r'))  # TODO: change the undump behavior here
    value_Matrix=pickle.load(file('pickleDump3.dump','r'))  # TODO: change the undump behavior here
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))  # TODO: change the undump behavior here
    cumulative=[]
    rei_UP=[]
    for SP in Uniprots:
        rei_UP.append(NodeID2MatrixNumber[SP])
    for i in range(0,iterations):
        random.shuffle(rei_UP)
        re_UP=rei_UP[:sample_size+1]
        source_MatrixID=re_UP[0]
        sinks_MatrixIDs=re_UP[0:]
        Currents=Info_circulation_for_Single_Node(conductance_Matrix,source_MatrixID,sinks_MatrixIDs,epsilon=1e-10)
        random.shuffle(rei_UP)
        targets=rei_UP[:sample_size/10]
        distances = dijkstra(value_Matrix, indices=source_MatrixID, unweighted=True)
        for index in targets:
            cumulative.append((distances[index],Currents[index,0]))
        for i in range(0,sample_size/5):
            index=np.argmin(distances)
            print distances
            cumulative.append((distances[index],Currents[index,0]))
            distances[index]=100
    pickle.dump(cumulative,file('Silverality.dump','w'))  # TODO: change the undump behavior here
    return cumulative


def analyze_Silverality():
    Silverality_List=pickle.load(file('Silverality.dump','r'))  # TODO: change the undump behavior here
    SuperDict={}
    for distance, value in Silverality_List:
        if distance not in SuperDict.keys():
            SuperDict[distance]=[]
        SuperDict[distance].append(value)
    
    for key, val in SuperDict.iteritems():
        print key, val
    
    ArrayDict={}
    StatsDict={}
    x=np.zeros((10,1))
    y1=np.zeros((10,1))
    y2=np.zeros((10,1))
    y3=np.zeros((10,1))
    i=0
    Srtd=sorted(SuperDict.iteritems(),key=operator.itemgetter(0))
    for key, val in Srtd:
        if i>9:
            break
        ArrayDict[key]=np.array(val)
        StatsDict[key]=(np.mean(ArrayDict[key]),np.std(ArrayDict[key]))
        x[i,0]=key
        y1[i,0]=StatsDict[key][0]
        y2[i,0]=StatsDict[key][0]+StatsDict[key][1]
        y3[i,0]=StatsDict[key][0]-StatsDict[key][1]
        i+=1

    
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


def Perform_Testing_Routines():
    raise NotImplementedError


def Compute_and_Store_circulation(handle, List_of_UPs, mongo_db_collection):
    print 'entering computation for: ', handle, 'with', len(List_of_UPs), 'UPs'
    NodeID2MatrixNumber, MatrixNumber2NodeID, ID2displayName, ID2Type, ID2Localization, Uniprots = pickle.load(file('pickleDump2.dump','r'))   # TODO: change the undump behavior here
    re_UP_List=[]
    for elt in List_of_UPs:
        if elt in NodeID2MatrixNumber.keys():
            re_UP_List.append(NodeID2MatrixNumber[elt])
    Informativity_Array=compute_sample_circulation_intensity_minimal(re_UP_List)   # TODO: change the undump behavior here
    UPs_Set=pickle.dumps(set(List_of_UPs))
    post={'GO_ID':handle,'UP_Set':UPs_Set,'Info_Array':pickle.dumps(Informativity_Array)}   # TODO: change the undump behavior here
    mongo_db_collection.insert(post)
    return Informativity_Array


def Compute_ponderated_info_circulation(UPs_2Binding_Affs):
    raise NotImplementedError

if __name__ == "__main__":
    characterise_eigenvects()


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


# Shut down HiNT analysis => Slightly improves the result

# Synchronious eigenvectors approach: protect agains entering into a forbidden list the target node
# start iterating matrix multiplications starting from the node1 to go to the node2
# enter each node visited in the forbidden set, except for node2
# terminate iterating when there are no more new reaches for node2 after all the interations

# Percentage of information reaching a given node compared to all the information reaching the node: eigenvalue approach too.
# Error we do: compute three times

# Ok, what is going on is that we have collections of ~ 300 elements completely screwing our system

# The problem that a information broadcasting between the elements of the same group is not a good thing, but a direct broadcasting into a reaction is actually
# what we need in our matrix.


# In order to be precise, we should not only take in account the power of bindinb between a molecule and protein and criticality of the protein, but also the abundance of the
# protein in the reactome

# => Done with the aboundance retrieval

# DONE: use sparse matrixes routines to calculate the number of connex elements in the graph
#   Problem: there are 58 disconnected sets.
#   Solution: retrieve the Node Ids of the main connex Set and write them into the neo4j graph, then retrieve only them

# DONE: markup of the major connex graph within neo4j database
#    Waiting for the execution


# DONE: calculate the distance graph
    # seems to work pretty well with Djikistra.
    # Can we perform a retrieval of specific nodes within distance X of the main component? 

# DONE: buid jump tables to compute the number of reactional transitions
#    Implemented by using djikstra algo from scipy.sparse.csgraph
#

# DONE: retrieve Pamela silver's degradation of the data with the time
#    Waiting for the execution
#    

# DONE: pull in the annotations regarding the proteins aboundances
#    
#    

# DONE: pull in the 300 essential targets from the EBI dude (John Overington)
#     Results aren't so conclusive. It seems that the protein concentration defenitely plays some role in the determining if a protein is a 
#     Target of an existing drug or not, butthe informativity seems not. Probably this is due to the fact that the targeted proteins are often 
#     cellular receptors.

# DONE: perform a localization factor pull-out for the Uniprots based on their proteins of attachement
#        Waiting for the execution

# DONE: broadcast to uniprots for the localization of the pointed proteins
