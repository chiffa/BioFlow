import os
from datetime import datetime

smtp_logging = False
smtp_logging_parameters = {
    'local_host': 'lpdpc4.epfl.ch',
    'local_mail_account': 'andrei@lpdpc4.epfl.ch',
    'reporting_target_mail': 'andrei.kucharavy@epfl.ch',
}


# TODO: add overrides from the command line by the user.
# TODO: reconcile with the servers.ini paths loading.
storage_location = os.path.join(os.environ['HOME'], 'bioflow')
sources_location = os.path.join(storage_location, 'sources')
output_location = os.path.join(storage_location, 'outputs')
internal_storage = os.path.join(storage_location, '.internal')
dumps_directory = os.path.join(internal_storage, 'dumps')
logs_directory= os.path.join(internal_storage, 'logs')


# build location used elsewhere
log_location = logs_directory
dump_location = dumps_directory

# build output directory
current_run_start_time = str(datetime.now())
current_run_start_time = current_run_start_time.replace(':', '.')
output_location = os.path.join(output_location, 'run started on ' + current_run_start_time)


# TODO: add overrides from the command line by the user.
# TODO: wrap it all in an "environment" class
env_skip_reactome = False
env_skip_hint = False
env_skip_biogrid = False
# CURRENTPASS: [BKI normalization] Make sure those are injected into the BioKnowledgeInterface
#  properly
env_use_background = True
env_bki_filter = ['biological_process']
env_bki_correlation_factors = (1, 1)
env_bki_ultraspec_clean = True  # BKI cleans ultra-specific terms
env_bki_ultraspec_lvl = 3  # How many proteins at most would be annotated by a term considered as
                           # ultra-specific


# those are deep configurations and should not be touched unless you know what you are doing:
switch_to_splu = False  # switching this to True incurs approximately an 100-fold slowdown
share_solver = True  # switching this to False incurs approximately a 50-fold slowdown
sparse_analysis_threshold = 200

# those are mostly debug flags and should not be touched
implicitely_threaded = True
psutil_main_loop_memory_tracing = False  # controls the log_mem behavior in conduction_routines.py
memory_source_allowed = False
node_current_in_debug = False

# Those are global variables that were used in order to check how algorithm variation performed
use_normalized_laplacian = False
fraction_edges_dropped_in_laplacian = 0.0

# This is the p-value that is set by user and is to be propagated throughout the analyses.
p_val_cutoff = 0.2
