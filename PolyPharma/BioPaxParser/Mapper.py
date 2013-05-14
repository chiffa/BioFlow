'''
Created on 12 mai 2013

@author: User
'''

from DBLoader import TableBuilder as TB
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

Reac2Unip_src = file('C:\Users\User\Documents\UCSD\Parsing_Reactome\uniprot_2_pathways.stid.txt','r')

i=0
while True:
    line = Reac2Unip_src.readline()
    i+=1
    words=line.split('\t')
    if len(words)>2:
        ACNUM=words[0]
        REACTId=words[1]
        session.add(TB.ReactomeId2acnum(ACNUM,REACTId))
    else:
        print i, line
    if not line:
        break

session.commit()


