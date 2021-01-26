import os
from datetime import datetime

smtp_logging = False
smtp_logging_parameters = {
    'local_host': 'lpdpc4.epfl.ch',
    'local_mail_account': 'andrei@lpdpc4.epfl.ch',
    'reporting_target_mail': 'andrei.kucharavy@epfl.ch',
}


# TODO: add overrides from the command line by the user.
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
skip_reactome = False
skip_hint = False
skip_biogrid = False

# those are deep configurations and should not be touched unless you know what you are doing:
switch_to_splu = False  # switching this to True incurs approximately an 100-fold slowdown
share_solver = True  # switching this to False incurs approximately a 50-fold slowdown
sparse_analysis_threshold = 200

# those are mostly debug flags and should not be touched
single_threaded = False
psutil_main_loop_memory_tracing = False  # controls the log_mem behavior in conduction_routines.py
memory_source_allowed = False
node_current_in_debug = False

# Those are global variables that were used in order to check how algorithm variation performed
use_normalized_laplacian = False
fraction_edges_dropped_in_laplacian = 0.0

# This is the p-value that is set by user and is to be propagated throughout the analyses.
p_val_cutoff = 0.05
