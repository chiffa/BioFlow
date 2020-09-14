import os

storage_location = os.path.join(os.environ['HOME'], 'bioflow')
sources_location = os.path.join(storage_location, 'sources')
output_location = os.path.join(storage_location, 'outputs')
internal_storage = os.path.join(storage_location, '.internal')
dumps_directory = os.path.join(internal_storage, 'dumps')
logs_directory = os.path.join(internal_storage, 'logs')

# TODO: add overrides from the command line by the user.

skip_reactome = False
skip_hint = False
skip_biogrid = False