"""
:created: 12 mai 2013
:@author: Andrei Kucharavy

Holds all the configurations of the environmental variables for the whole project
"""
__author__ = 'ank'

import pickle
from pprint import PrettyPrinter
from os import path, makedirs
from pymongo import MongoClient
from BioFlow.Utils.ConfigsIO import parse_configs, conf_file_path_flattener
from BioFlow.configs.internals_config import edge_type_filters, Leg_ID_Filter, fudge, Adjacency_Martix_Dict, Conductance_Matrix_Dict

Servers, Options, Sources, Predictions = parse_configs()

# SERVERS CONFIGURATIONS
MongoDB_url = Servers['PRODUCTION']['mongodb_server']
neo4j_server = Servers['PRODUCTION']['server_neo4j']
ReadSourceDBs = conf_file_path_flattener(Sources)

# REQUIRED PARAMETERS
GeneOntology = ReadSourceDBs['GO']
ReactomeBioPax = ReadSourceDBs['REACTOME']
UNIPROT_source = ReadSourceDBs['UNIPROT']
Hint_csv = ReadSourceDBs['HINT']  # attention, for me it is tab-separated
BioGRID = ReadSourceDBs['BIOGRID']

# OPTIONAL PARAMETERS
Chromosome_source = ReadSourceDBs['CHROMOSOMES']
Chromosome_file_filter = Sources['CHROMOSOMES']['namepattern']
# PROTEIN ABOUNDANCES?

# Defines Mongodb properties and connections
pymongo_prefix = Sources['INTERNAL']['mongoprefix']
pymongo_suffix = Sources['INTERNAL']['mongosuffix']

# Builds static entry points for mongodatabase access
client = MongoClient(MongoDB_url)
db = client.PolyPharma_database
UP_rand_samp = db[pymongo_prefix+"UP_r_samples"+pymongo_suffix]
Interactome_rand_samp = db[pymongo_prefix+"Interactome_samples"+pymongo_suffix]


#########################################################################
#  Defines the dumping locations for different intermediate computations
#########################################################################
class Dumps(object):
    """
    A class that contains and controls all the dumps related to accelerated loading of mappings between the graph DB
    and the mapping matrix holders
    """
    prefix = str(path.abspath(path.dirname(__file__)+'/dumps'))
    prefix_2 = Sources['INTERNAL']['dumpprefix']
    postfix = '.dump'

    if not path.isdir(prefix+prefix_2):
        makedirs(prefix+prefix_2)

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

    RNA_seq_counts_compare = prefix + prefix_2 + '/RNA_seq_compare'+postfix


####################################################
#  Defines the locations to output actual results  #
####################################################
class Outputs(object):
    prefix = str(path.abspath(path.dirname(__file__)+'/outputs'))
    GO_GDF_output = prefix + '/GO_Analysis_output.gdf'
    Interactome_GDF_output = prefix + '/Interactome_Analysis_output.gdf'
    RNA_pre_filter_output = prefix + '/RNA_pre_filter_output.tsv'
    cross_refs = prefix + '/cross_refs.tsv'


#########################################################################
#  Declares overloaded IDs, pickles from the dumps of already computed  #
#########################################################################
IDFilter = []
if path.isfile(Dumps.Forbidden_IDs):
    IDFilter = pickle.load(file(Dumps.Forbidden_IDs, 'r'))

#####################################################################################
# Where the RNA counts source, hits and background deduced from it are to be found  #
#####################################################################################
# TODO: these should be input dynamically
RNA_source = "/home/ank/Documents/External_Predictions/Ben_RNA_seq/counts.tsv"
Hits_source = "/home/ank/projects_files/2015/Hung_Ji_essential_genes/shortlist.csv"
Background_source = "/home/ank/projects_files/2014/Poly_Pharma/HJ-screen/Allgene_R2.csv"

# these are remappings to the inner master database IDs
prename1 = Hits_source[:-4]+'_'+'pPh_name_maps.txt'
prename2 = Hits_source[:-4]+'_'+'pPh_id_list.csv'
bgList = Background_source[:-4]+'_'+'pPh_id_list.csv'

if __name__ == "__main__":
    pp = PrettyPrinter(indent=4)
    pp.pprint((Servers, Options, Sources))
    pass