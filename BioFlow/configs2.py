"""
Holds all the configurations of the environmental variables for the whole project
"""
import os
import pickle
from os import path, makedirs
from pprint import PrettyPrinter

from pymongo import MongoClient

from BioFlow.utils.ConfigsIO import parse_configs, conf_file_path_flattener
from BioFlow.utils.general_utils import high_level_os_io as SF

Servers, Options, Sources, Predictions = parse_configs()

# SERVERS CONFIGURATIONS
mongo_db_url = Servers['PRODUCTION']['mongodb_server']
neo4j_server = Servers['PRODUCTION']['server_neo4j']
ReadSourceDBs = conf_file_path_flattener(Sources)

# REQUIRED PARAMETERS
GeneOntology = ReadSourceDBs['GO']
ReactomeBioPax = ReadSourceDBs['REACTOME']
UNIPROT_source = ReadSourceDBs['UNIPROT']
Hint_csv = ReadSourceDBs['HINT']  # attention, for me it is tab-separated
BioGRID = ReadSourceDBs['BIOGRID']

# TAXONOMY IDs for Uniprot parsing:
up_tax_ids = [tax_id.strip() for tax_id in Sources['UNIPROT']['tax_ids'].split(',') if
              tax_id not in ('', ' ')]

# OPTIONAL PARAMETERS
Chromosome_source = ReadSourceDBs['CHROMOSOMES']
Chromosome_file_filter = Sources['CHROMOSOMES']['namepattern']
# PROTEIN ABOUNDANCES?

# Defines Mongodb properties and connections
pymongo_prefix = Sources['INTERNAL']['mongoprefix']
pymongo_suffix = Sources['INTERNAL']['mongosuffix']

# Builds static entry points for mongodatabase access
client = MongoClient(mongo_db_url)
db = client.PolyPharma_database
UP_rand_samp = db[pymongo_prefix+"UP_r_samples"+pymongo_suffix]
Interactome_rand_samp = db[pymongo_prefix+"Interactome_samples"+pymongo_suffix]


class Dumps(object):
    """
    A class that contains and controls all the dumps related to accelerated loading of mappings
    between the graph DB and the mapping matrix holders
    """
    prefix = path.join(path.abspath(
        path.join(path.dirname(__file__), os.pardir)), 'dumps')
    SF.mkdir_recursive(prefix)

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


class Outputs(object):
    """
    Defines the locations to output actual results
    """
    prefix = path.join(path.abspath(
        path.join(path.dirname(__file__), os.pardir)), 'outputs')
    SF.mkdir_recursive(prefix)

    GO_GDF_output = prefix + '/GO_Analysis_output.gdf'
    Interactome_GDF_output = prefix + '/Interactome_Analysis_output.gdf'
    RNA_pre_filter_output = prefix + '/RNA_pre_filter_output.tsv'
    cross_refs = prefix + '/cross_refs.tsv'


#  Declares overloaded IDs, pickles from the dumps of already computed
IDFilter = []
if path.isfile(Dumps.Forbidden_IDs):
    IDFilter = pickle.load(file(Dumps.Forbidden_IDs, 'r'))


# Where the RNA counts source, hits and background deduced from it are to be found  #
# TODO: these should be input dynamically; or at least from a different source file because of a
# different modification frequency
RNA_source = "/home/ank/Documents/External_Predictions/Ben_RNA_seq/counts.tsv"
Hits_source = "/home/andrei/support/tmp/Chr_10.txt"
Background_source = "/home/ank/projects_files/2014/Poly_Pharma/HJ-screen/Allgene_R2.csv"

# these are remappings to the inner master database IDs
prename1 = Hits_source[:-4] + '_' + 'pPh_name_maps.txt'
prename2 = Hits_source[:-4] + '_' + 'pPh_id_list.csv'
bgList = Background_source[:-4] + '_' + 'pPh_id_list.csv'


output_location = path.join(path.abspath(
    path.join(path.dirname(__file__), os.pardir)), 'outputs')


if __name__ == "__main__":
    pp = PrettyPrinter(indent=4)
    pp.pprint((Servers, Options, Sources))
