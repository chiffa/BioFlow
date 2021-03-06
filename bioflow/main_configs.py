"""
Holds all the configurations of the environmental variables for the whole project
"""
import os
import pickle
from os import path, makedirs
from pprint import PrettyPrinter
from bioflow.configs_manager import parse_config, compute_full_paths
from bioflow.user_configs import output_location, dump_location, log_location
from bioflow.utils.general_utils import high_level_os_io as hl_os_io
from datetime import datetime


current_run_start_time = str(datetime.now())
current_run_start_time = current_run_start_time.replace(':', '.')


# # TODO: remap those to the values in user configs
# dump_location = dumps_directory
# output_location = path.join(path.abspath(os.path.expanduser('~')), 'outputs' + '_' +
#                             current_run_start_time)
# # TODO: [run path refactor] Move this to auto_analyze (TOP)
# output_location = path.join(output_location, 'run started on ' + current_run_start_time)
# # TODO: [run path refactor] Move this to auto_analyze (TOP)
# log_location = logs_directory  # TODO: this does nothing, given that logs already imported

# TODO: [run path refactor] Move this to auto_analyze (TOP)
hl_os_io.mkdir_recursive(dump_location)
hl_os_io.mkdir_recursive(output_location)
hl_os_io.mkdir_recursive(log_location)

Servers = parse_config('servers')
Sources = parse_config('sources')
DB_locations = parse_config('online_dbs')

# SERVERS CONFIGURATIONS
mongo_db_url = os.getenv('MONGOURL', Servers['PRODUCTION']['mongodb_server'])
neo4j_server = os.getenv('NEO4URL', Servers['PRODUCTION']['server_neo4j'])
source_db_paths = compute_full_paths(Sources, DB_locations, Servers['PRODUCTION'])

# Locations of data source files to be passed to parsers
gene_ontology_path = source_db_paths['GO']
reactome_biopax_path = source_db_paths['REACTOME']
uniprot_path = source_db_paths['UNIPROT']
hint_csv_path = source_db_paths['HINT']  # attention, for me it is tab-separated
biogrid_path = source_db_paths['BIOGRID']

organism_meta_flag = Sources['META']['organism']

if organism_meta_flag == 'Human':
    # marbach_path = source_db_paths['MARBACH']
    # marbach_mode = Sources['MARBACH']['mode']
    # marbach_sig = Sources['MARBACH']['significance']
    trrust_path = source_db_paths['TRRUST']
    trrust_sig = Sources['TRRUST']['significance']
    # cellnet_path = source_db_paths['CELLNET']
    # cellnet_sig = Sources['CELLNET']['significance']
    complexes_path = source_db_paths['COMPLEXPORTAL']

else:
    trrust_path = None
    trrust_sig = None
    complexes_path = None

phosphosite_path = source_db_paths['PHOSPHOSITE']
phosphosite_organism = Sources['PHOSPHOSITE']['organism']

# TAXONOMY IDs for Uniprot parsing:
up_tax_ids = [tax_id.strip() for tax_id in Sources['UNIPROT'][
    'tax_ids'].split(',') if tax_id not in ('', ' ')]

# Defines Mongodb properties and connections
pymongo_prefix = Sources['INTERNAL']['mongoprefix']
pymongo_suffix = Sources['INTERNAL']['mongosuffix']
estimated_comp_ops = int(Sources['INTERNAL']['compops'])  # pairwise flows computed per second


class Dumps(object):
    """
    A class that contains and controls all the dumps related to accelerated loading of mappings
    between the graph DB and the mapping matrix holders
    """
    prefix = dump_location

    prefix_2 = Sources['INTERNAL']['dumpprefix']
    postfix = '.dump'

    if not path.isdir(prefix + prefix_2):
        makedirs(prefix + prefix_2)

    ###################################
    matrix_LS = prefix + prefix_2 + '/dump5' + postfix
    interactome_maps = prefix + prefix_2 + '/dump2' + postfix
    eigen_VaMat = prefix + prefix_2 + '/eigen_valmat.csv'
    eigen_ConMat = prefix + prefix_2 + '/eigen_conmat.csv'
    val_eigen = prefix + prefix_2 + '/pickleDump' + postfix
    cond_eigen = prefix + prefix_2 + '/pickleDump_0_5' + postfix
    interactome_adjacency_matrix = prefix + prefix_2 + '/pickleDump3' + postfix
    interactome_laplacian_matrix = prefix + prefix_2 + '/pickleDump4' + postfix
    UniP_att = prefix + prefix_2 + '/UP_Attach' + postfix
    Main_Connex_group = prefix + prefix_2 + '/Connex_group' + postfix
    Forbidden_IDs = prefix + prefix_2 + '/ForbiddenIDs' + postfix
    Adj_degree = prefix + prefix_2 + '/Adjacency_degree.csv'
    Silverality = prefix + prefix_2 + '/Silverality' + postfix
    InfoArray = prefix + prefix_2 + '/sample_array' + postfix
    Interactome_Analysis_memoized = prefix + \
        prefix_2 + '/Interactome_memoization' + postfix

    Up_dict_dump = prefix + prefix_2 + '/Uniprot_dict' + postfix
    GO_dump = prefix + prefix_2 + '/GO' + postfix
    GO_builder_stat = prefix + prefix_2 + '/GO_builder_stats' + postfix
    GO_Mats = prefix + prefix_2 + '/GO_mats' + postfix
    GO_Infos = prefix + prefix_2 + '/GO_Infos' + postfix
    GDF_debug = prefix + prefix_2 + '/GDF_debug.gdf'
    GO_Inflated = prefix + prefix_2 + '/GO_inflated' + postfix
    GO_Analysis_memoized = prefix + prefix_2 + '/GO_memoization' + postfix
    GO_Indep_Linset = prefix + prefix_2 + '/GO_Indep_linset' + postfix

    RNA_seq_counts_compare = prefix + prefix_2 + '/RNA_seq_compare' + postfix

    # TODO: this has nothing to do with in the configs either
    analysis_set_display_names = prefix + prefix_2 + '/current_analysis_set_name_maps.txt'
    analysis_set_bulbs_ids = prefix + prefix_2 + '/current_analysis_set_bulbs_id_list.csv'
    background_set_bulbs_ids = prefix + prefix_2 + '/current_background_set_bulbs_id_list.csv'


# TODO: modifier insertion in this case is pretty effed up I would guess - the evaluation is done
#  immediately to build the variables
# TODO: [run path refactor] pipe hdd save destination here (TOP)
class Outputs(object):
    """
    Defines the locations to output actual results
    """
    prefix = output_location

    GO_GDF_output = prefix + '/GO_Analysis_output.gdf'
    Interactome_GDF_output = prefix + '/Interactome_Analysis_output.gdf'
    RNA_pre_filter_output = prefix + '/RNA_pre_filter_output.tsv'
    cross_refs = prefix + '/cross_refs.tsv'

    knowledge_network_stats = prefix + '/knowledge_network_stats.png'
    interactome_network_stats = prefix + '/interactome_network_stats.png'

    knowledge_network_output = prefix + '/knowledge_stats.tsv'
    interactome_network_output = prefix + '/interactome_stats.tsv'


# CURRENTPASS:
class NewOutputs(object):

    def __init__(self, modifier=''):

        if modifier != '':
            root_path = path.join(output_location, modifier)
        else:
            root_path = output_location

        hl_os_io.mkdir_recursive(root_path)

        self.GO_GDF_output = path.join(root_path, 'GO_Analysis_output.gdf')
        self.Interactome_GDF_output = path.join(root_path, 'Interactome_Analysis_output.gdf')
        self.RNA_pre_filter_output = path.join(root_path, 'RNA_pre_filter_output.tsv')
        self.cross_refs = path.join(root_path, 'cross_refs.tsv')

        self.knowledge_network_stats = path.join(root_path, 'knowledge_network_stats.png')
        self.interactome_network_stats = path.join(root_path, 'interactome_network_stats.png')

        self.knowledge_network_output = path.join(root_path, 'knowledge_analysis_stats.tsv')
        self.interactome_network_output = path.join(root_path, 'interactome_analysis_stats.tsv')

        self.interactome_network_scatterplot = path.join(root_path, 'interactome.png')


#  Declares overloaded IDs, pickles from the dumps of already computed
forbidden_neo4j_ids = []
if path.isfile(Dumps.Forbidden_IDs):
    try:
        forbidden_neo4j_ids = pickle.load(open(Dumps.Forbidden_IDs, 'rb'))
    except Exception as e:
        print('exception encountered: %s' % (e))

# TODO: cut this away
# Where the RNA counts bioflow, hits and background deduced from it are to be found  #
# these are defaults that can be overriden by changing parameters to "cast analysis set" function
#  from neo4j db io module
rna_source = "/home/ank/Documents/External_Predictions/Ben_RNA_seq/counts.tsv"
analysis_protein_ids_csv = "/home/andrei/support/tmp/Chr_10.txt"
background_protein_ids_csv = "/home/ank/projects_files/2014/Poly_Pharma/HJ-screen/Allgene_R2.csv"

if __name__ == "__main__":
    pp = PrettyPrinter(indent=4)
    pp.pprint((Servers, Sources))
