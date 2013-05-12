'''
Created on 12 mai 2013

@author: User
'''

from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()
