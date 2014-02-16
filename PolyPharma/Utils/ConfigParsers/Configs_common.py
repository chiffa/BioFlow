__author__ = 'ank'

import os

db_root='/scratch/local_home/'
source_ext_db_roots='/home/ank/Documents/External_DBs_Store/'
source_ext_prediction_roots='/home/ank/Documents/External_Predictions/'

rootdir = os.path.abspath(os.path.join(os.path.dirname(__file__),'../../configs/'))
shortnames = ['servers','options','sources','predictions']
configsfiles = [rootdir + '/' + name + '.ini' for name in shortnames ]