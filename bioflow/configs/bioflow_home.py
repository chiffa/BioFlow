"""
Core start-up loop that is required by all the other modules to operate.
Non-logged, so the configuration is minimal
"""
import os
from datetime import datetime

# figure where to read configs from and store the logs
bioflow_home_directory = os.path.join(os.environ['HOME'], 'bioflow')

if "BIOFLOWHOME" in os.environ:
    bioflow_home_directory = os.environ["BIOFLOWHOME"]


sources_location = os.path.join(bioflow_home_directory, 'sources')
output_location = os.path.join(bioflow_home_directory, 'outputs')
internal_storage = os.path.join(bioflow_home_directory, '.internal')
dumps_directory = os.path.join(internal_storage, 'dumps')
logs_directory = os.path.join(internal_storage, 'logs')
confs_location = os.path.join(bioflow_home_directory, 'configs')

# build location used elsewhere
log_location = logs_directory
dump_location = dumps_directory

# build output directory
current_run_start_time = str(datetime.now())
current_run_start_time = current_run_start_time.replace(':', '.')
output_location = os.path.join(output_location, 'run started on ' + current_run_start_time)
