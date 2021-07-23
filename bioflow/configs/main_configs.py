"""
Holds all the configurations of the environmental variables for the whole project
"""
import os
import pickle
import yaml
from os import path, makedirs
from pprint import pprint
from collections import defaultdict

from bioflow.configs.configs_manager import compute_full_paths
from bioflow.configs.bioflow_home import output_location, log_location, dump_location, \
    sources_location, confs_location, internal_storage
from bioflow.configs.internal_configs import reactome_forbidden_nodes, uniprot_forbidden_nodes
from bioflow.utils.general_utils import high_level_os_io as hl_os_io
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


# create all the relevant folders
hl_os_io.mkdir_recursive(dump_location)
hl_os_io.mkdir_recursive(output_location)
hl_os_io.mkdir_recursive(log_location)
hl_os_io.mkdir_recursive(confs_location)


# # potentially ensure that the base config.yaml is findable upon pip installation
# configs_rootdir = os.path.abspath(
#         os.path.dirname(__file__))
# log.info(configs_rootdir)


# deploy configs if don't exist
user_yaml_confs = os.path.join(confs_location, 'main_configs.yaml')
hl_os_io.copy_if_doesnt_exist('configs.yaml', user_yaml_confs)


# read the configs that the user might have defined
log.info('loading configs from %s' % user_yaml_confs)
with open(user_yaml_confs, 'r', encoding='utf8') as user_yaml_configs:
    configs_loaded = yaml.safe_load(user_yaml_configs)
    Servers = configs_loaded['Servers']
    Sources = configs_loaded['Sources']
    # REFACTOR: [better confs]: potential improvement - parse all organisms and select here
    DB_locations = configs_loaded['DB_locations']
    organism = Sources['META']['organism']
    user_settings = configs_loaded['User_Settings']
    log.info('loaded configs from %s' % organism)


# SERVERS CONFIGURATIONS
mongo_db_url = os.getenv('MONGOURL', Servers['mongodb_server'])
neo4j_server_url = os.getenv('NEO4URL', Servers['neo4j_server'])
neo4j_user = 'neo4j'
neo4j_user = os.getenv('NEO4USER', Servers['neo4j_user'])
source_db_paths = compute_full_paths(Sources, DB_locations, sources_location)


# Locations of data source files to be passed to parsers
gene_ontology_path = source_db_paths['GO']
reactome_biopax_path = source_db_paths['REACTOME']
uniprot_path = source_db_paths['UNIPROT']
hint_csv_path = source_db_paths['HINT']  # attention, for me it is tab-separated
biogrid_path = source_db_paths['BIOGRID']

if organism == 'Human':
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
neo4j_db_name = Sources['INTERNAL'].get('neo4jdb', None)


class Dumps(object):
    """
    A class that contains and controls all the dumps related to accelerated loading of mappings
    between the graph DB and the mapping matrix holders
    """
    prefix = os.path.join(dump_location, organism)
    postfix = '.dump'

    if not path.isdir(prefix):
        makedirs(prefix)

    ###################################
    matrix_LS = os.path.join(prefix, 'dump5.dump')
    interactome_maps = os.path.join(prefix, 'dump2.dump')
    eigen_VaMat = os.path.join(prefix, 'eigen_valmat.csv')
    eigen_ConMat = os.path.join(prefix, 'eigen_conmat.csv')
    val_eigen = os.path.join(prefix, 'pickleDump.dump')
    cond_eigen = os.path.join(prefix, 'pickleDump_0_5.dump')
    interactome_adjacency_matrix = os.path.join(prefix, 'pickleDump3.dump')
    interactome_laplacian_matrix = os.path.join(prefix, 'pickleDump4.dump')
    UniP_att = os.path.join(prefix, 'UP_Attach.dump')
    Main_Connex_group = os.path.join(prefix, 'Connex_group.dump')
    Forbidden_IDs = os.path.join(prefix, 'ForbiddenIDs.dump')
    Adj_degree = os.path.join(prefix, 'Adjacency_degree.csv')
    Silverality = os.path.join(prefix, 'Silverality.dump')
    InfoArray = os.path.join(prefix, 'sample_array.dump')
    Interactome_Analysis_memoized = os.path.join(prefix, 'Interactome_memoization.dump')

    Up_dict_dump = os.path.join(prefix, 'Uniprot_dict.dump')
    GO_dump = os.path.join(prefix, 'GO.dump')
    GO_builder_stat = os.path.join(prefix, 'GO_builder_stats.dump')
    GO_Mats = os.path.join(prefix, 'GO_mats.dump')
    GO_Infos = os.path.join(prefix, 'GO_Infos.dump')
    GDF_debug = os.path.join(prefix, 'GDF_debug.gdf')
    GO_Inflated = os.path.join(prefix, 'GO_inflated.dump')
    GO_Analysis_memoized = os.path.join(prefix, 'GO_memoization.dump')
    GO_Indep_Linset = os.path.join(prefix, 'GO_Indep_linset.dump')

    RNA_seq_counts_compare = os.path.join(prefix, 'RNA_seq_compare.dump')

    # those are temporary storage of cast sets of IDs and backgrounds
    analysis_set_display_names = prefix + '/current_analysis_set_name_maps.txt'
    analysis_set_bulbs_ids = prefix + '/current_analysis_set_bulbs_id_list.csv'
    background_set_bulbs_ids = prefix + '/current_background_set_bulbs_id_list.csv'


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
        self.knowledge_network_scatterplot = path.join(root_path, 'knowledge.png')

        # cluster_complement

        self.knowledge_clusters_output = path.join(root_path, 'knowledge_clusters_stats.tsv')
        self.interactome_clusters_output = path.join(root_path, 'interactome_clusters_stats.tsv')

        self.interactome_clusters_scatterplot = path.join(root_path, 'interactome_clusters.png')
        self.knowledge_clusters_scatterplot = path.join(root_path, 'knowledge_clusters.png')


# This can now be reweighted by policy modification
laplacian_default_type_edge_weighting = {
    "is_part_of_collection": 0.5,
    "is_same": 100,
    "is_catalysant": 1,
    "is_reaction_participant": 1,
    "is_part_of_complex": 1,
    "is_regulant": 1,
    "is_interacting": 1,
    "is_weakly_interacting": 0.5,
    "is_likely_same": 1,
    "is_able_to_modify": 1,
}

laplacian_default_source_edge_weighting = defaultdict(lambda: 1)
# we don't care about sources by default


adjacency_default_type_edge_weighting = {
    "is_part_of_collection": 0.5,
    "is_same": 1,
    "is_catalysant": 0.33,
    "is_reaction_participant": 0.33,
    "is_part_of_complex": 0.33,
    "is_regulant": 0.33,
    "is_interacting": 0.33,
    "is_weakly_interacting": 0.15,
    "is_likely_same": 0.1,
    "is_able_to_modify": 0.1,
}

adjacecency_default_source_edge_weighting = defaultdict(lambda: 1)
# we don't care about sources by default


#  Declares overloaded IDs, pickles from the dumps of already computed
forbidden_neo4j_ids = []
if path.isfile(Dumps.Forbidden_IDs):
    try:
        forbidden_neo4j_ids = pickle.load(open(Dumps.Forbidden_IDs, 'rb'))
    except Exception as e:
        log.critical('exception encountered: %s' % (e))
        raise e


# REFACTOR: [better confs]: wrap the 'env_' variables it all in an "environment" wrap
smtp_logging = bool(user_settings['smtp_logging'])
smtp_logging_parameters = user_settings['smtp_logging_parameters']

env_skip_reactome = bool(user_settings['environment']['skip_reactome'])
env_skip_hint = bool(user_settings['environment']['skip_hint'])
env_skip_biogrid = bool(user_settings['environment']['skip_biogrid'])

env_use_background = bool(user_settings['environment']['use_background'])
env_bki_filter = user_settings['environment']['bki_filter']
env_bki_correlation_factors = tuple(user_settings['environment']['bki_correlation_factors'])
env_bki_ultraspec_clean = bool(user_settings['environment']['bki_ultraspec_clean'])
# BKI cleans ultra-specific terms
env_bki_ultraspec_lvl = int(user_settings['environment']['bki_ultraspec_lvl'])
# How many proteins at most would be annotated by a term considered as

switch_to_splu = bool(user_settings['solver']['switch_to_splu'])
# switching this to True incurs approximately an  100-fold slowdown
share_solver = bool(user_settings['solver']['share_solver'])
# switching this to False incurs approximately a 50-fold slowdown
line_loss = float(user_settings['solver']['line_loss'])
# This is the line loss for the approximate matrix inversion - basically the fudge for cholesky

implicitely_threaded = bool(user_settings['debug_flags']['implicitely_threaded'])
psutil_main_loop_memory_tracing = bool(user_settings['debug_flags']['psutil_main_loop_memory_tracing'])
# controls the log_mem behavior in conduction_routines.py
memory_source_allowed = bool(user_settings['debug_flags']['memory_source_allowed'])
node_current_in_debug = bool(user_settings['debug_flags']['node_current_in_debug'])

use_normalized_laplacian = user_settings['use_normalized_laplacian']
fraction_edges_dropped_in_laplacian = user_settings['fraction_edges_dropped_in_laplacian']

sparse_analysis_threshold = int(user_settings['analysis']['sparse_analysis_threshold'])
default_background_samples = int(user_settings['analysis']['default_background_samples'])
default_p_val_cutoff = float(user_settings['analysis']['default_p_val_cutoff'])
min_nodes_for_p_val = int(user_settings['analysis']['min_nodes_for_p_val'])

neo4j_autobatch_threshold = int(configs_loaded['Servers'].get('neo4j_autobatch_threshold', 5000))


if env_skip_hint:
    laplacian_default_source_edge_weighting['HINT'] = 0
    adjacecency_default_source_edge_weighting['HINT'] = 0

if env_skip_biogrid:
    laplacian_default_source_edge_weighting['BioGRID'] = 0
    adjacecency_default_source_edge_weighting['BioGRID'] = 0

if env_skip_reactome:
    laplacian_default_source_edge_weighting['Reactome'] = 0
    adjacecency_default_source_edge_weighting['Reactome'] = 0


if __name__ == "__main__":
    pass
    # pprint((Servers, Sources))


