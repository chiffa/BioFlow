'''
Created on Feb 25, 2013

@author: akucahravy
'''


from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker



Base = declarative_base()

lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb', echo=False)

Base.metadata.reflect(lite_engine)
meta=Base.metadata

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

dropTable='--'
#dropTable='UNIPROT2GO'
#dropTable='UNIPROTid2acnum'
#dropTable='UNIPROT_Prot'

for tbl in reversed(meta.sorted_tables):
    print tbl
    if tbl.key == dropTable:
        print "dropping",dropTable
        tbl.drop(lite_engine)


