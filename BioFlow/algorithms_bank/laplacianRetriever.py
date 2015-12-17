import hashlib
import json
import pickle
import random
import string
from collections import defaultdict
from copy import copy
from csv import reader
from itertools import combinations, chain
from math import log
from pprint import PrettyPrinter
from random import shuffle
from time import time
from warnings import warn

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import shortest_path

from BioFlow.algorithms_bank import conduction_routines as CR
from BioFlow.main_configs import Dumps, Outputs, UP_rand_samp, background_set_bulbs_ids
from BioFlow.molecular_network.InteractomeInterface import MatrixGetter
from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph
from BioFlow.utils.gdfExportInterface import GdfExportInterface
from BioFlow.utils.io_Routines import dump_object, undump_object


class LaplacianRetriever(object):

    def __init__(self):
        pass
