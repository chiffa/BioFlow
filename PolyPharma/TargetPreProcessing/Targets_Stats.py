'''
Created on Apr 5, 2013

@author: akucahravy
'''

from DBLoader import TableBuilder as TB
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

EffectiveDrugLIGs = ['AZM','CPF','4AX','KAN','NMY','PAR','RPT','RFP','RBT','SRY','SAN','TOY','VAN']

def get_eff_drugs_Targets():
    query=session.query(TB.DRUG2UNIPROT).\
                    filter(TB.DRUG2UNIPROT.drugLig.in_(EffectiveDrugLIGs)).\
                    all()
    
    retdict={}
    for row in query:
        if str(row.uniprotid) not in retdict.keys():
            retdict[str(row.uniprotid)]=0
        retdict[str(row.uniprotid)]+=1
            
    return retdict