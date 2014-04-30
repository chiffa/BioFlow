"""
:created: 12 mai 2013
:@author: Andrei Kucharavy

Holds all the configurations of the environmental variables for the whole project
"""
__author__='ank'

from PolyPharma.Utils.ConfigParsers.Configs_parser import parse_configs, sourcefile_compilator
from pprint import PrettyPrinter
from pymongo import MongoClient
from os import path
import pickle
import os.path


Servers, Options, Sources, Predictions = parse_configs()


SQLite_location = Servers['PRODUCTION']['local_sqlite']
MongoDB_url = Servers['PRODUCTION']['mongodb_server']
neo4j_server = Servers['PRODUCTION']['server_neo4j']
ReadSourceDBs = sourcefile_compilator(Sources)


GeneOntology = ReadSourceDBs['GO']
ReactomeBioPax = ReadSourceDBs['REACTOME']
UNIPROT_source = ReadSourceDBs['UNIPROT']
Hint_csv = ReadSourceDBs['HINT']  #attention, for me it is tab-separated
Protein_aboundances = ReadSourceDBs['ABOUNDANCES']
sedEffFileName = ReadSourceDBs['SIDER']  #TODO: improve mappings from Drugs to secondary effects.
Chromosome_source = ReadSourceDBs['CHROMOSOMES']
Chromosome_file_filter = Sources['CHROMOSOMES']['namepattern']

ReadSourcePredictions = sourcefile_compilator(Predictions)


# Targets file is assumed to be tab-separated, with the first column containing fuzzy versions of "names" of "gene names"
# from uniport and the next three - information relative to binding to this target.
Targets_File = ReadSourcePredictions['NEFLANAVIR']
# File from which to load the names of the 300 most frequent targets
Targets_File2 = ReadSourcePredictions['OVERINGTON']


# ExactDict is the dictionnary used to perform a precise matching between the fuzzy target names and the SwissProt IDs required
# for a lookup in the database
from TargetPreProcessing.neflanavir_parser import subdict
Targets_dict = subdict
from TargetPreProcessing.Overington_parser import subdict2
Targets_dict2 = subdict2

################################################
#  Defines MongeDb properties and connections
################################################
# pymongo_prefix = "human_"
# pymongo_prefix = "mice_"
pymongo_prefix = "yeast_"
pymongo_suffix = "_v_2"

client = MongoClient(MongoDB_url)
db = client.PolyPharma_database
UP_rand_samp = db[pymongo_prefix+"UP_r_samples"+pymongo_suffix]
Interactome_rand_samp = db[pymongo_prefix+"Interactome_samples"+pymongo_suffix]


#######################################################################
#  Defines how much confidence we have into the different interactions
#######################################################################
# Refers to the groups of links between the nodes that should be treated in the same manner
edge_type_filters = {
    "Group" : ["is_part_of_collection"],                                  # Group relation group
    "Same" : ["is_same"],                                                 # Same relation group
    "Reaction" : ["is_Catalysant", "is_reaction_particpant"],             # Reaction relation group
    "Contact_interaction" : ["is_part_of_complex", "is_Regulant"],        # Contact_interaction relation group
    "HiNT_Contact_interaction" : ["is_interacting"],                      # Contact_interaction relation group
    "possibly_same" : ["is_possibly_same"],
    }


# Coefficients values for the value_Matrix
Adjacency_Martix_Dict = {"Group":0.5,
             "Same":1,
             "Reaction":0.33,
             "Contact_interaction":0.33,
             "possibly_same":0.1,
             }


# Coefficients values for the conductance_Matrix
Conductance_Matrix_Dict = {"Group":0.5,
             "Same":100,
             "Reaction":1,
             "Contact_interaction":1,
             "possibly_same":0.1,
             }

#########################################################################
#  Defines the dumping locations for different intermediate computations
#########################################################################
class Dumps(object):
    """
    A class that contains and controls all the dumps related to accelerated loading of mappings between the graph DB
    and the mapping matrix holders
    """
    prefix = str(path.abspath(path.dirname(__file__)+'/dumps'))
    # TODO: achtung:explosive here
    prefix_2 = '/yeast'
    postfix = '.dump'

    if not os.path.isdir(prefix+prefix_2):
        os.makedirs(prefix+prefix_2)

    ###################################
    matrix_LS = prefix + prefix_2 + '/dump5'+postfix
    matrix_corrs = prefix + prefix_2 + '/dump2'+postfix
    eigen_VaMat = prefix + prefix_2 + '/eigen_valmat.csv'
    eigen_ConMat = prefix + prefix_2 + '/eigen_conmat.csv'
    val_eigen = prefix + prefix_2 + '/pickleDump'+postfix
    cond_eigen = prefix + prefix_2 + '/pickleDump_0_5'+postfix
    ValMat = prefix + prefix_2 + '/pickleDump3'+postfix
    ConMat = prefix + prefix_2 + '/pickleDump4'+postfix
    UniP_att = prefix + prefix_2 + '/UP_Attach'+postfix
    Main_Connex_group = prefix + prefix_2 + '/Connex_group'+postfix
    Forbidden_IDs = prefix + prefix_2 + '/ForbiddenIDs'+postfix
    Adj_degree = prefix + prefix_2 + '/Adjacency_degree.csv'
    Silverality = prefix + prefix_2 + '/Silverality'+postfix
    InfoArray = prefix + prefix_2 + '/sample_array'+postfix
    Interactome_Analysis_memoized = prefix + prefix_2 + '/Interactome_memoization'+postfix

    Up_dict_dump = prefix + prefix_2 + '/Uniprot_dict'+postfix
    GO_dump = prefix + prefix_2 + '/GO'+postfix
    GO_builder_stat = prefix + prefix_2 + '/GO_builder_stats'+postfix
    GO_Mats = prefix + prefix_2 + '/GO_mats'+postfix
    GO_Infos = prefix + prefix_2 + '/GO_Infos'+postfix
    GDF_debug = prefix + prefix_2 + '/GDF_debug.gdf'
    GO_Inflated = prefix + prefix_2 + '/GO_inflated'+postfix
    GO_Analysis_memoized = prefix + prefix_2 + '/GO_memoization'+postfix
    GO_Indep_Linset = prefix + prefix_2 + '/GO_Indep_linset'+postfix


#########################################################################
#  Defines the locations to output actual results
#########################################################################
class Outputs(object):
    prefix = str(path.abspath(path.dirname(__file__)+'/outputs'))
    GO_GDF_output = prefix + '/GO_Analysis_output.gdf'
    Interactome_GDF_output = prefix + '/Interactome_Analysis_output.gdf'


###########################################################################################
#  Defines what nodes are to be masked to avoid conduction ovrload of non-informative nodes
###########################################################################################
Leg_ID_Filter = ['H+', 'ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2', 'NTP',]

IDFilter = []
if os.path.isfile(Dumps.Forbidden_IDs):
    IDFilter = pickle.load(file(Dumps.Forbidden_IDs,'r'))

##########################################################################
#  Fundge for matrix diagolizations of matrixes and other solver functions
##########################################################################
fudge = 1e-10

if __name__ == "__main__":
    pp=PrettyPrinter(indent=4)
    pp.pprint((Servers,Options,Sources))
