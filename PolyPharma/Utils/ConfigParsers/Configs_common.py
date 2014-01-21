__author__ = 'ank'

import os

rootdir = os.path.abspath(os.path.join(os.path.dirname(__file__),'../../configs/'))
shortnames = ['servers','options']
configsfiles = [rootdir + '/' + name + '.ini' for name in shortnames ]