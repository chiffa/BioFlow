'''
Created on Feb 25, 2013

@author: akucahravy
'''
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy import Index
from sqlalchemy import ForeignKey

lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb',echo=False)

Base = declarative_base()

# TODO: bring GORapidfire Table in this declaration sheet
# TODO: pull it out from the current package and into a highest-level file,
# just as file name references

class GO_Term(Base):
    '''
    Encodes all the informations relative to the GO terms and their
    definitions
    '''
    __tablename__ = 'GO_Term'
    
    termid = Column(String, primary_key=True)
    name = Column(String)
    namespace = Column(String)
    definition = Column(String)
    # Normally theere should also be an additional field for the term informativity
    # But it will be added in a separate table
    # The columns for the computation-intensive tables should be indexed too
    # The caching should be done much more agressively, before any communication with database
    
    
    
    def __init__(self,termid,name,namespace,definition):
        self.termid = termid
        self.name = name
        self.namespace = namespace
        self.definition = definition
    
    def __repr__(self):
        return "<GO_Term('%s','%s','%s','%s')>" % (self.termid,
                                                   self.name,
                                                   self.namespace,
                                                   self.definition)

class GO_Informativity(Base):
    '''
    Encodes the GO terms informativity.
    Normally should be added as an additional column to the 
    GO_Term table, but was forgotten and later added as a separate table
    In order not to break the incoming references from the 
    '''
    __tablename__='GO_Informativity'
    
    go=Column(String,ForeignKey('GO_Term.termid'),primary_key=True,)
    informativity=Column(Float)
    
    def __init__(self, go, informativity):
        self.go=go
        self.informativity=informativity
    
    def __repr__(self):
        return "<GO_Informativity('%s','%s')>"%(self.go,
                                           self.informativity)


class GO_Relation(Base):
    '''
    Encodes the relations between the GO terms
    '''
    __tablename__ = 'GO_Relation'
    
    relid = Column(Integer, primary_key=True)
    gofrom = Column(String,ForeignKey('GO_Term.termid'))
    reltype = Column(String)
    goto = Column(String,ForeignKey('GO_Term.termid'))
    
    def __init__(self, gofrom, reltype, goto):
        self.gofrom = gofrom
        self.reltype = reltype
        self.goto = goto
    
    def __repr__(self):
        return "<GO_Relation('%s','%s','%s')>" % (self.gofrom,
                                                       self.reltype,
                                                       self.goto)

class UNIPROT_Prot(Base):
    '''
    Encodes the informations regarding a specific protein;
    uses it's SWISSPROT ID as a unique identifier
    '''
    __tablename__ = 'UNIPROT_Prot'
    
    uniprotid = Column(String, primary_key=True)
    full_name = Column(String)
    orf_name = Column(String)
    taxid = Column(String)
    #taxid seems to be unique for each uniprot protein
    
    def __init__(self,uniprotid,full_name,orf_name,taxid):
        self.uniprotid = uniprotid
        self.full_name = full_name
        self.orf_name = orf_name
        self.taxid = taxid
    
    def __repr__(self):
        return "<UNIPROT('%s','%s','%s','%s')>" % (self.uniprotid,
                                                   self.full_name,
                                                   self.orf_name,
                                                   self.taxid)
    
class GeneName2UNIPROT(Base):
    '''
    Encodes the relation between the gene names and the UNIPROT IDs
    '''
    __tablename__='GENENAME2UNIPROTID'
    
    gene2uniprotid=Column(Integer,primary_key=True)
    genename=Column(String)
    uniprotid=Column(String,ForeignKey('UNIPROT_Prot.uniprotid'))
    
    def __init__(self,genename,uniprotid):
        self.genename=genename
        self.uniprotid=uniprotid
        
    def __repr__(self):
        return "<Gene2Uniprot('%s','%s')>" % (self.genename,
                                              self.uniprotid)


class PDB2UNIPROT(Base):
    '''
    Do not try to join it with anything else!
    There are items in uniprot id which are not in UNIPROT_PROT.uniprotID!!!!
    '''
    
    __tablename__='PDB2UNIPROTID'
    
    gene2uniprotid=Column(Integer,primary_key=True)
    pdbid=Column(String, index=True)
    uniprotid=Column(String, ForeignKey('UNIPROT_Prot.uniprotid'), index=True,)
    
    def __init__(self,pdbid,uniprotid):
        self.pdbid=pdbid
        self.uniprotid=uniprotid
        
    def __repr__(self):
        return "<Pdb2Uniprot('%s','%s')>" % (self.pdbid,
                                             self.uniprotid)

   
class Drug(Base):
    '''
    Encodes the information regarding the drugs used in the 
    MYCTU dataset, based on their Ligand short name
    # TODO: fuse it with drugname table and make the reference depend
    # On the drug Name, rather then the lig Code. Manage the LIG
    # Code redundancies by multiplying the number of the LIG-columns
    # Alternative: retrieve LIG codes for all the drugs...
    '''
    __tablename__='DRUG'
    
    LIG=Column(String,primary_key=True)
    name=Column(String)
    
    def __init__(self,LIG,name):
        self.LIG=LIG
        self.name=name
        
    def __repr__(self):
        return "<DRUG('%s','%s')>" % (self.LIG,
                                      self.name)
        

class DrugName(Base):
    '''
    Table referencing all the drugs taken in account while working on the
    Human drug-target interaction dataset
    '''
    __tablename__='DRUGNAME'
    
    name=Column(String,primary_key=True)
    
    def __init__(self,name):
        self.name=name
        
    def __repr__(self):
        return "<DRUGName('%s')>" % (self.name)


class DRUG2UNIPROT(Base):
    '''
    Relation between drugs and uniprot identifiers extracted from
    the Pr. Bourne's work dataset of Mycobacterium Tuberculosis
    '''
    __tablename__ = 'DRUG2UNIPROT'
    
    updrugrelid = Column(Integer, primary_key=True)
    uniprotid = Column(String,ForeignKey('UNIPROT_Prot.uniprotid'))
    drugLig = Column(String,ForeignKey('GO_Term.termid'))
    energy = Column(Float)
    
    def __init__(self,uniprotid,drugLig,energy):
        self.uniprotid = uniprotid
        self.drugLig = drugLig
        self.energy = energy
    
    def __repr__(self):
        return "<DRUG2UNIPROT('%s','%s','%s')>" % (
                                                   self.uniprotid,
                                                   self.drugLig,
                                                   self.energy)

class MedDRASecEffTerms(Base):
    '''
    Table encoding the MedDRA secondary effects term names
    '''
    __tablename__='MedDRATerms'
    
    TermID=Column(String,primary_key=True)
    name=Column(String)
    
    def __init__(self,TermID,name):
        self.TermID=TermID
        self.name=name
        
    def __repr__(self):
        return "<MedDRATerms('%s','%s')>" % (self.TermID,
                                             self.name)

class DRUG2HUMUNIPROT(Base):
    '''
    Relation between drugs and uniprot identifiers extracted from
    the Pr. Bourne's work dataset on human secondary effects
    '''
    __tablename__ = 'DRUG2HUMUNIPROT'
    
    updrugrelid2 = Column(Integer, primary_key=True)
    uniprotid = Column(String,ForeignKey('UNIPROT_Prot.uniprotid'))
    drugName = Column(String,ForeignKey('DRUGNAME.name'))
    energy = Column(Float)
    
    def __init__(self,uniprotid,drugName,energy):
        self.uniprotid = uniprotid
        self.drugName = drugName
        self.energy = energy
    
    def __repr__(self):
        return "<DRUG2HUMUNIPROT('%s','%s','%s')>" % (
                                                   self.uniprotid,
                                                   self.drugName,
                                                   self.energy)
 
class DRUG2SECEFF(Base):
    '''
    relation between druygs and secondary effects, as inferred from teh database 
    downloaded from sideffects.embl.de
    '''
    __tablename__ = 'DRUG2SECEFF'
    
    drugsecrelid = Column(Integer, primary_key=True)
    drugName = Column(String,ForeignKey('DRUGNAME.name'))
    SecEffTerm = Column(String,ForeignKey('MedDRATerms.TermID'))
    
    def __init__(self,drugName,SecEffTerm):
        self.drugName = drugName
        self.SecEffTerm = SecEffTerm
    
    def __repr__(self):
        return "<DRUG2SECEFF('%s','%s')>" % (
                                             self.drugName,
                                             self.SecEffTerm)
    
class UNIPROT2GO(Base):
    '''
    relation between the uniprot terms and 
    '''
    __tablename__ = 'UNIPROT2GO'
    
    upgorelid = Column(Integer, primary_key=True)
    uniprotid = Column(String,ForeignKey('UNIPROT_Prot.uniprotid'))
    goid = Column(String,ForeignKey('GO_Term.termid'))
    
    def __init__(self,uniprotid,goid):
        self.uniprotid = uniprotid
        self.goid = goid
    
    def __repr__(self):
        return "<UNIPROT2GO('%s','%s')>" % (
                                              self.uniprotid,
                                              self.goid)

class UNIPROTid2acnum(Base):
    __tablename__ = 'UNIPROTid2acnum'
    
    upacrelid = Column(Integer, primary_key=True)
    uniprotid = Column(String,ForeignKey('UNIPROT_Prot.uniprotid'))
    acnum = Column(String)
    
    def __init__(self,uniprotid,acnum):
        self.uniprotid = uniprotid
        self.acnum = acnum
    
    def __repr__(self):
        return "<UNIPROTid2acnum('%s','%s')>" % (
                                              self.uniprotid,
                                              self.acnum)

class UNIPROT_FullGO2(Base):
    __tablename__ = 'UNIPROT_FullGO2'
    
    UPFGOrelid2 = Column(Integer, primary_key=True)
    UNIPROTID = Column(String, ForeignKey('UNIPROT_Prot.uniprotid'),index=True)
    timeout= Column(Integer,index=True)
    GOs = Column(String, ForeignKey('GO_Term.termid'),index=True)
    
    Index('CombinedUPIDtimeout',UNIPROTID,timeout)
    
    def __init__(self, UNIPROTID, timeout, GOs):
        self.UNIPROTID = UNIPROTID
        self.timeout = timeout
        self.GOs = GOs
    
    def __repr__(self):
        return "<UNIPROT_FullGO2('%s','%s','%s')>" % (self.UNIPROTID,
                                                      self.timeout,
                                                      self.GOs)

class UNIPROT_SpecInf(Base):
    __tablename__ = 'UNIPROT_SpecInf'
    
    UPFSpecInfrelid = Column(Integer, primary_key=True)
    UNIPROTID = Column(String, ForeignKey('UNIPROT_Prot.uniprotid'),index=True)
    timeout= Column(Integer,index=True)
    TotalInfo = Column(Float,index=True)
    TotalGOs = Column(Integer,index=True)
    
    Index('CombinedUPFSpecInfrel1',UNIPROTID,timeout)
    Index('CombinedUPFSpecInfrel2',UNIPROTID,timeout,TotalInfo)
    Index('CombinedUPFSpecInfrel3',UNIPROTID,timeout,TotalGOs)
    
    def __init__(self, UNIPROTID, timeout, TotalInfo, TotalGOs):
        self.UNIPROTID = UNIPROTID
        self.timeout = timeout
        self.TotalInfo = TotalInfo
        self.TotalGOs = TotalGOs
    
    def __repr__(self):
        return "<UNIPROT_SpecInf('%s','%s','%s','%s')>" % (self.UNIPROTID,
                                                           self.timeout,
                                                           self.TotalInfo,
                                                           self.TotalGOs)


'''
<============================================>
'''
        
Base.metadata.create_all(lite_engine)
