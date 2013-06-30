'''
Created on Mar 19, 2013
@author: akucahravy
'''

from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import operator
import math
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


class Protein:
    '''
    A wrapper object for Proteins, allowing to conviniently store information about
    the proteins itself and it's relation with the other GO terms
    
    Has fields:
            name:            the name of the protein
            importance:      relative importance of the protein in an external dataset
            UpStream:        *Objects* (GOs, Source, Sink) that are *upstream* in the information flow
            DownStream:      *Objects* (GOs, Source, Sink) that are *downstream* in the information flow
            DownStreamInf:   Sum of Informativities of the objects that are DownStream
            InfConduct:      Information Conductivity of the protein, terms behind it and proteins behind the terms
            Flow:            Flow of information through this protein
            Finalize:        Checks if the InfConduct have been computed for this GOTerm
    
    '''
    name=''
    importance=-1
    UpStream=[]
    DownStream=[]
    DownStreamInf=-1
    InfConduct=-1
    Flow=0
    Finalized=False
    isLeft=True
    
    def __init__(self,name,importance,UpStream=[],DownStream=[]):
        '''
        Class Constructor:
        
        Arguments:  name = the UNIPROT_ID of the protein
                    importance = the external importance score for the protein
                    left = List of objects (GOterm) to which the object
                            is connected to it's left
                    right = List of objects (GOterm) to which the object
                            is connected to it's right
                            
        Attention! The left and right arguments should be lists of *Objects*, not *Strings*
                    The way to go is to first create the object with void left and right
                    And then fill them with add_to_right and add_to_left functions
        
        Attention! Three class attributes are not added in the beginning because they
                    rely heavily on the structure of the graph and have to be calculated
                    either with the "Finalize()" method or in the GOFlowSystem
        '''
        self.name=name
        self.importance=importance
        self.UpStream=UpStream[:]
        self.DownStream=DownStream[:]
    
    def add_to_UpStream(self,Object):
        '''
        Adds the object to the UpStream
        '''
        self.UpStream.append(Object)
    
    def add_to_DownStream(self,Object):
        '''
        Adds the object to the DownStream
        '''
        self.DownStream.append(Object)
    
    def ShowChildren(self,Object):
        for obj in self.DownStream:
            print obj
    
    def Finalize(self):
        '''
        This method should be called once the structure setting down is finished
        and just before proceeding to the computation of the information flow in 
        the GoFlowSystem container
        '''
        if self.DownStream==[]:
            self.InfConduct=self.importance # limit behavior if conductivity of source
            self.Finalized=True             # is infinite
        else:
            totalInformativity=0
            for Obj in self.DownStream:
                if not Obj.Finalized:
                    Obj.Finalize()
                totalInformativity+=Obj.InfConduct
            self.DownStreamInf=totalInformativity
            self.InfConduct=self.importance
            self.Finalized=True
    
    def __repr__(self):
        return "<Protein('%s','%s','%s','%s','%s','%s','%s')>" % (self.name,
                                                                  self.importance,
                                                                  self.InfConduct,
                                                                  self.Flow,
                                                                  self.Finalized,
                                                                  self.isLeft,
                                                                  len(self.DownStream))

class GOTerm:
    '''
    A wrapper object for GOTerms, allowing to conviniently store information about
    them and their relations with the Proteins  itself and it's relation with the other GO terms
    
    Has fields:
            TermID:          the name of the protein
            informativity:   the informativity of the GO term
            UpStream:        *Objects* (Proteins) that are *upstream* in the information flow
            DownStream:      *Objects* (Proteins) that are *downstream* in the information flow
            DownStreamInf:   Sum of Informativities of the objects that are DownStream
            InfConduct:      Information Conductivity of the GO term and proteins behind it
            Flow:            Flow of information through this GOTerm
            Finalize:        Checks if the InfConduct have been computed for a this GOTerm
    
    '''
    TermID=''
    informativity=-1
    UpStream=[]
    DownStream=[]
    DownStreamInf=-1
    InfConduct=-1
    Flow=0
    Finalized=False
    
    def __init__(self, TermID, informativity, UpStream=[], DownStream=[]):
        '''
        Class Constructor:
        
        Arguments:  TermID = the UNIPROT_ID of the protein
                    informativity = the external importance score for the protein
                    left = List of objects (Proteins) to which the object
                            is connected to it's left
                    right = List of objects (Proteins) to which the object
                            is connected to it's right
                            
        Attention! The left and right arguments should be lists of *Objects*, not *Strings*
                    The way to go is to first create the object with void left and right
                    And then fill them with add_to_right and add_to_left functions
        
        Attention! Three class attributes are not added in the beginning because they
                    rely heavily on the structure of the graph and have to be calculated
                    either with the "Finalize()" method or in the GOFlowSystem
        '''
        self.TermID=TermID
        self.informativity=informativity
        self.UpStream=UpStream[:]
        self.DownStream=DownStream[:]
    
    def add_to_UpStream(self,Object):
        '''
        Adds the object to the UpStream
        '''
        self.UpStream.append(Object)
    
    def add_to_DownStream(self,Object):
        '''
        Adds the object to the DownStream
        '''
        self.DownStream.append(Object)
    
    def ShowChildren(self,Object):
        for obj in self.DownStream:
            print obj
        
    def Finalize(self):
        '''
        This method should be called once the structure setting down is finished
        and just before proceeding to the computation of the information flow in 
        the GoFlowSystem container
        '''
        if self.DownStream==[]:
            self.InfConduct=self.informativity # limit behavior if conductivity of source
            self.Finalized=True             # is infinite. Used to compute the conductivity
                                            # In case of single-sided similarity computation 
        else:
            totalInformativity=0
            for Obj in self.DownStream:
                if not Obj.Finalized:
                    Obj.Finalize()
                totalInformativity+=Obj.InfConduct
            self.InfConduct=1.0/(1.0/float(self.informativity)+1.0/float(totalInformativity))
            self.Finalized=True              

    def __repr__(self):
        return "<GOTerm('%s','%s','%s','%s','%s','%s')>" % (self.TermID,
                                                            self.informativity,
                                                            self.InfConduct,
                                                            self.Flow,
                                                            self.Finalized,
                                                            len(self.DownStream))

    ###
    # TODO: Check what GOs are common between the GO lists
    # form the objects for the GOs and the proteins
    
    ###
    # TODO: Check that all the proteins without GO annotation have
    #been eliminated from the set

class GOFlowSystem:
    '''
    A Wrapper functions used to calculate information flow through a set of GOs
    General Idea:
    
                |  ProteinsLeft  |  GOTerms  |  ProteinsRight  |
    
                                |---->GO1--->|
                |--->Protein2-->|            |
    Source ---->|               |-|-->GO2--->|----->Protein3------->Sink
                |--->Protein1---->|          |
                                  |-->GO3--->|
    '''
    # Lists where the key is the name (string) and the value is the object
    ObjProtLeft={} 
    ObjGOs={}
    ObjProtRight={}
    GO_Type='null'
    timeout=3
    EdgeList={}
    corrfactor=math.exp(1)

    
    def __init__(self,ProteinsLeft,ProteinsRight={},GO_Type='null',timeout=3,corrfactor=math.exp(1)):
        '''
        Initializes all the relations between the protein and the GO terms
        given an initial protein list and properly filled database from the GOAnnot package
        
        Arguments:  ProteinsRight = Python dictionary of the type {'Uniprot_ID' : Confidence_lvl}

                    ProteinsLeft = Python dictionary of the type {'Uniprot_ID' : Confidence_lvl}
                              if == {} (default), the Left Side is replaced by a supersink
                                       (protein with no importance, but annotated with all GOs )

                    GO_Type = 'null' by default
                         or = 'biological_process','cellular_component','molecular_function'
        
                    timeout = 3 by default. 
                              The timeout on the "regulates" transition from the
                                GO term database
        '''
        self.GO_Type=GO_Type
        self.timeout=timeout
        
        tempGOs1=[]
        tempGOs2=[]
        
        for ProteinName in ProteinsLeft.keys():
            ProtObj=Protein(ProteinName,ProteinsLeft[ProteinName])
            ProtObj.isLeft=True
            self.ObjProtLeft[ProteinName]=ProtObj
            tempGOs1=tempGOs1+self.queryGO(ProteinName)
        
        if ProteinsRight != {}:
            #The most complex case
            for ProteinName in ProteinsRight.keys():
                ProtObj=Protein(ProteinName,ProteinsRight[ProteinName])
                ProtObj.isLeft=False
                self.ObjProtRight[ProteinName]=ProtObj
                tempGOs2=tempGOs2+self.queryGO(ProteinName)
        
            finalGOs=list(set(tempGOs1) & set(tempGOs2))
        
            for GO in finalGOs:
                self.ObjGOs[GO]=GOTerm(GO,self.getInfoGO(GO))
            
        else: 
            for GO in tempGOs1:
                self.ObjGOs[GO]=GOTerm(GO,self.getInfoGO(GO))
            
    def buildTopology(self):
        
        # For the sake of method stability we have to eliminate all the proteins without
        # annotation within the set of the GO terms we consider as "common"
        
        UnannotatedProts=[]
        
        for ProteinName in self.ObjProtLeft.keys():
            GOs=self.queryGO(ProteinName)
            Eliminated=True
            for GO in GOs:
                if GO in self.ObjGOs.keys():
                    Eliminated=False
                    self.ObjProtLeft[ProteinName].add_to_DownStream(self.ObjGOs[GO])
                    self.ObjGOs[GO].add_to_UpStream(self.ObjProtLeft[ProteinName])
            if Eliminated:
                UnannotatedProts.append(ProteinName)
        
        for ProtName in UnannotatedProts:
            del self.ObjProtLeft[ProtName]
        
        if len(UnannotatedProts)>0:
            print 'deleted unannotated left proteins:', UnannotatedProts
        
        UnannotatedProts=[]
        
        for ProteinName in self.ObjProtRight.keys():
            GOs=self.queryGO(ProteinName)
            Eliminated=True
            for GO in GOs:
                if GO in self.ObjGOs.keys():
                    Eliminated=False
                    self.ObjProtRight[ProteinName].add_to_UpStream(self.ObjGOs[GO])
                    self.ObjGOs[GO].add_to_DownStream(self.ObjProtRight[ProteinName])
            if Eliminated:
                UnannotatedProts.append(ProteinName)
        
        for ProtName in UnannotatedProts:
            del self.ObjProtRight[ProtName]
            
        if len(UnannotatedProts)>0:
            print 'deleted unannotated right proteins:', UnannotatedProts
        
        # Finalize Everything: Normally the last line should suffice, but just in case...
        
        for ProteinName in self.ObjProtRight.keys():
            self.ObjProtRight[ProteinName].Finalize()
        for GOName in self.ObjGOs.keys():
            self.ObjGOs[GOName].Finalize()
        for ProteinName in self.ObjProtLeft.keys():
            self.ObjProtLeft[ProteinName].Finalize()
            
    def ComputeFlow(self,TotalTension):
        # Now we Finally can compute the information flow through the terms, starting from
        # the left and finishing on the right
        
        # correlate the Tension with the total importance of proteins, so that the protein
        # similarity score can be expressed as a the value of the flow
        
        # transform GO annotation values so that they specifically mediate the 
        # informativity of a given protein 
        
        for ProteinObj in self.ObjProtLeft.values():
            ProteinObj.Flow=TotalTension*ProteinObj.InfConduct
            GOTotInf=0
            # Let's split the information flow through a protein through the DownStream GO terms
            for GO in ProteinObj.DownStream:
                GOTotInf+=GO.InfConduct
            GOTotInf=float(GOTotInf)
            for GO in ProteinObj.DownStream:
                GO.Flow+=ProteinObj.Flow/GOTotInf*GO.InfConduct
                self.EdgeList[(ProteinObj.name,GO.TermID)]=ProteinObj.Flow/GOTotInf*GO.InfConduct
        
        # So, now we've got the information conductivity for each of the LeftProteins and 
        # GO Terms. Let's do the RightTerms\
        
        for GO in self.ObjGOs.values():
            ProtTotInf=0.0
            for ProteinObj in GO.DownStream:
                ProtTotInf+=ProteinObj.InfConduct
            ProtTotInf=float(ProtTotInf)
            for ProteinObj in GO.DownStream:
                ProteinObj.Flow+=GO.Flow/ProtTotInf*ProteinObj.InfConduct
                self.EdgeList[(GO.TermID,ProteinObj.name)]=GO.Flow/ProtTotInf*ProteinObj.InfConduct
       
    def queryGO(self,ProteinName):
        
        Query=[]
        if self.GO_Type not in ['biological_process','cellular_component','molecular_function']:
            Query= session.query(TableBuilder.UNIPROT_FullGO2).\
                            filter(TableBuilder.UNIPROT_FullGO2.UNIPROTID==ProteinName).\
                            filter(TableBuilder.UNIPROT_FullGO2.timeout==self.timeout).\
                            all()
        else:
            Query= session.query(TableBuilder.UNIPROT_FullGO2).\
                            join(TableBuilder.GO_Term).\
                            filter(TableBuilder.UNIPROT_FullGO2.UNIPROTID==ProteinName).\
                            filter(TableBuilder.UNIPROT_FullGO2.timeout==self.timeout).\
                            filter(TableBuilder.GO_Term.namespace==self.GO_Type).\
                            all()
        
        result = []
        
        for element in Query:
            result.append(str(element.GOs))
        
        result=list(set(result))   
        return result

    def getInfoGO(self,GO):
        query=session.query(TableBuilder.GO_Informativity).\
                        filter(TableBuilder.GO_Informativity.go == GO).\
                        one()
        return float( math.log(float(query.informativity)+0.0001)/math.log(self.corrfactor))
    
    def Verify(self):
        LeftCurrent=0.0
        GoCurrent=0.0
        RightCurrent=0.0
        for prot in self.ObjProtRight.values():
            RightCurrent+=prot.Flow
        for prot in self.ObjProtLeft.values():
            LeftCurrent+=prot.Flow
        for GoObj in self.ObjGOs.values():
            GoCurrent+=GoObj.Flow
        
        return [LeftCurrent,GoCurrent,RightCurrent]
    
    def sortTargets(self):
        '''
        Returns the final processed list of Targets sorted according to importance
        And the GO list of annotation under the format: 
        '''
        LeftProt2Val={}
        ID2GoVal={}
        RightProt2Val={}
        
        for ProtName in self.ObjProtLeft.keys():
            LeftProt2Val[ProtName]=self.ObjProtLeft[ProtName].Flow
        
        for ProtName in self.ObjProtRight.keys():
            RightProt2Val[ProtName]=self.ObjProtRight[ProtName].Flow        
        
        for GoId in self.ObjGOs.keys():
            ID2GoVal[GoId]=self.ObjGOs[GoId].Flow
        
        sortedLP = sorted(LeftProt2Val.iteritems(),key=operator.itemgetter(1),reverse=True)
        sortedRP = sorted(RightProt2Val.iteritems(),key=operator.itemgetter(1),reverse=True)
        sortedGO = sorted(ID2GoVal.iteritems(),key=operator.itemgetter(1),reverse=True)
        
        return [sortedLP,sortedGO,sortedRP]

    def getDetailedInformationGO(self,Go):
        TermName = str( session.query(TableBuilder.GO_Term).\
                                filter(TableBuilder.GO_Term.termid==Go).\
                                one().name )
        LinkingIn={elt:value for (elt,value) in self.EdgeList.items() if elt[1]==Go}
        LinkingOut={elt:value for (elt,value) in self.EdgeList.items() if elt[0]==Go}
        
        return [TermName,LinkingIn,LinkingOut]

    def getDetailedInformationProt(self,Prot,left):
        TermName = str( session.query(TableBuilder.UNIPROT_Prot).\
                                filter(TableBuilder.UNIPROT_Prot.uniprotid==Prot).\
                                one().full_name )
        LinkingIn={}
        LinkingOut={}
        if not left:
            LinkingIn={elt:value for (elt,value) in self.EdgeList.items() if elt[1]==Prot}
        if left:
            LinkingOut={elt:value for (elt,value) in self.EdgeList.items() if elt[0]==Prot}
        # TODO: add normalization regarding total GO power
        TotalGOPower=str(session.query(TableBuilder.UNIPROT_SpecInf).\
                                    filter(TableBuilder.UNIPROT_SpecInf.UNIPROTID==Prot).\
                                    filter(TableBuilder.UNIPROT_SpecInf.timeout==self.timeout).\
                                    one().TotalInfo)
        
        return [TermName,LinkingIn,LinkingOut,TotalGOPower]        
    
    def getDetailedReportGOs(self,firstN='all'):
        SortedGOs=self.sortTargets()[1]
        
        reportBuffer=''
        
        if firstN!='all' and firstN<len(SortedGOs):
            SortedGOs=SortedGOs[:firstN]

        GoCurrent=0.0
        for GoObj in self.ObjGOs.values():
            GoCurrent+=GoObj.Flow
        
        for GO in SortedGOs:
            reportBuffer+='GOGOGOGOGO>> '+str(GO[0])+' | '+str("{0:.2f}".format(self.ObjGOs[GO[0]].informativity))+' | '+str("{0:.2f}".format(GO[1]/GoCurrent*100))+'%\n'
            GOReport=self.getDetailedInformationGO(GO[0])
            name=GOReport[0]
            InEdges=GOReport[1]
            OutEdges=GOReport[2]
            reportBuffer+=name
            InSum=sum(InEdges.values())
            reportBuffer+='\n'+'Proteins from left  : '
            for (key,val) in sorted(InEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                reportBuffer+=str(key[0])+':'+str("{0:.2f}".format(val/InSum*100))+'% | '
            OutSum=sum(OutEdges.values())
            reportBuffer+='\n'+'Proteins from right : '
            for (key,val) in sorted(OutEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                reportBuffer+=str(key[1])+':'+str("{0:.2f}".format(val/OutSum*100))+'% | '
            reportBuffer+='\n'+'<<<<GOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGO<<'+'\n'
            
        return reportBuffer
            
    
    # TODO: separate reports for the right and left sides
    
    
    def getDetailedReportProts(self,firstNL='all',firstNR='all'):
        ST=self.sortTargets()
        LP=ST[0]
        RP=ST[2]
        
        reportBuffer=''
        
        #Processing LP, which is never empty:
        
        if firstNL!='all' and firstNL<len(LP):
            LP=LP[:firstNL]
        
        if firstNR!='all' and firstNR<len(RP):
            RP=RP[:firstNR]
        
        reportBuffer+='Processing left: \n ==================================== \n'
        
        for Prot in LP:
            reportBuffer+='PPPPPPP>> '+str(Prot[0])+' | '+str("{0:.2f}".format(self.ObjProtLeft[Prot[0]].DownStreamInf))+' | '+str("{0:.2f}".format(self.ObjProtLeft[Prot[0]].InfConduct))+' | '+str("{0:.2f}".format(Prot[1]))+'\n'
            ProtReport=self.getDetailedInformationProt(Prot[0],True)
            name=ProtReport[0]
            InEdges=ProtReport[1]
            OutEdges=ProtReport[2]
            MaxGOInf=ProtReport[3]
            reportBuffer+=name+' | '+MaxGOInf
            InSum=sum(InEdges.values())
            if InSum!=0:
                reportBuffer+='\n'+'GOs UpStream   : '
                for (key,val) in sorted(InEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                    reportBuffer+=str(key[0])+':'+str("{0:.2f}".format(val/InSum*100))+'% | '
            OutSum=sum(OutEdges.values())
            if OutSum!=0:
                reportBuffer+='\n'+'GOs DownStream : '
                for (key,val) in sorted(OutEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                    reportBuffer+=str(key[1])+':'+str("{0:.2f}".format(val/OutSum*100))+'% | '
            reportBuffer+='\n'+'<<<<<PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP<<'+'\n'
            
        reportBuffer+='Processing right \n ==================================== \n'
         
        for Prot in RP:
            reportBuffer+='>>>>>>>>> '+str(Prot[0])+' | '+str("{0:.2f}".format(self.ObjProtRight[Prot[0]].DownStreamInf))+' | '+str("{0:.2f}".format(self.ObjProtRight[Prot[0]].InfConduct))+' | '+str("{0:.2f}".format(Prot[1]))+'\n'
            ProtReport=self.getDetailedInformationProt(Prot[0],False)
            name=ProtReport[0]
            InEdges=ProtReport[1]
            OutEdges=ProtReport[2]
            MaxGOInf=ProtReport[3]
            reportBuffer+=name+' | '+MaxGOInf
            InSum=sum(InEdges.values())
            if InSum!=0:
                reportBuffer+='\n'+'GOs UpStream   : '
                for (key,val) in sorted(InEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                    reportBuffer+=str(key[0])+':'+str("{0:.2f}".format(val/InSum*100))+'% | '
            OutSum=sum(OutEdges.values())
            if OutSum!=0:
                reportBuffer+='\n'+'GOs DownStream : '
                for (key,val) in sorted(OutEdges.iteritems(),key=operator.itemgetter(1),reverse=True):
                    reportBuffer+=str(key[1])+':'+str("{0:.2f}".format(val/OutSum*100))+'% | '
            reportBuffer+='\n'+'<<<<PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP<<'+'\n'
         
            
        return reportBuffer
    
        # TODO: introduce the elimination of empty dictionaries, treating left and right
        # side separately to allow the self-similarity instructions

        
#<====================================================================>

###
# TODO: normalize the total informativity of proteins so it is the same value (say 1.0 )on the both sides of the 
# GOFlow system
# Idea: renormalize the importance of proteins so that their DownStreamInf are identical
# when they are tested in the "similarity" mode. i.g. normalize regarding all the GOs
# summed informativity for this protein and not only the the commonGOs with the other
# proteins

###
# TODO: perform a protein-wise normalization of the GO terms so that the proteins that 
# are highly annotated don't pull too much attention on them => it's more of a Zipf distribution question here

amp=3.330

dict1={'ZZEF1_HUMAN':amp,'ZZZ3_HUMAN':amp,'ZXDC_HUMAN':amp}
dict2={'ZXDB_HUMAN':amp,'ZWINT_HUMAN':amp}
dict3={}
dict4={'ZXDB_HUMAN':amp}

GFS=GOFlowSystem(dict1,dict1)

### 
# We can balance the distribution of importance between the GO annotation and protein
# Importance by increasing the proteins importance

###
# I really have an impression that the GOs should be renormalized, so that a molecular function tree
# has a zero informativity at the top, in the same way as the biological_process informativity

##
#
# percentage-wise channeling

GFS.buildTopology()
GFS.ComputeFlow(1000)
#print GFS.Verify()
#GFS.showSortedTargets()
#print GFS.getDetailedInformationGO('0008270')
#print GFS.getDetailedInformationProt('ZXDB_HUMAN')
#print GFS.getDetailedInformationProt('ZWINT_HUMAN')
print GFS.getDetailedReportGOs()
print GFS.getDetailedReportProts()