'''
Created on Mar 1, 2013

@author: akucahravy

!!!!!This sections should be entirely rebuild and merged with the DB builder in order
!!!!!to account for the "Foreing Key" domain

'''

from sqlalchemy import create_engine
from DBLoader.TableBuilder import GO_Term
from DBLoader.TableBuilder import UNIPROT_Prot
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy import ForeignKey
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Base = declarative_base()

class GO_Rapidfire(Base):
    __tablename__ = 'GO_Rapidfire'
    
    rapfireid = Column(Integer, primary_key=True)
    rootGO = Column(String)
    timeout= Column(Integer)
    otherGO = Column(String)
    
    def __init__(self, rootGO, timeout, otherGO):
        self.rootGO = rootGO
        self.timeout = timeout
        self.otherGO = otherGO
    
    def __repr__(self):
        return "<GO_Rapidfire('%s','%s','%s')>" % (self.rootGO,
                                                   self.timeout,
                                                   self.otherGO)
        
#class UNIPROT_FullGO(Base):
#    __tablename__ = 'UNIPROT_FullGO'
#    
#    UPFGOrelid = Column(Integer, primary_key=True)
#    UNIPROTID = Column(String) #Need a Foreign key, but the foreign key have to be placed in
#                                # The same file as the rest of the schema definition
#    timeout= Column(Integer)
#    GOs = Column(String) # Same thing here.
#    
#    def __init__(self, UNIPROTID, timeout, GOs):
#        self.UNIPROTID = UNIPROTID
#        self.timeout = timeout
#        self.GOs = GOs
#    
#    def __repr__(self):
#        return "<UNIPROT_FullGO('%s','%s','%s')>" % (self.UNIPROTID,
#                                                     self.timeout,
#                                                     self.GOs)

Base.metadata.create_all(lite_engine)