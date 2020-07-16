import os

storage_location = os.path.join(os.environ['HOME'], 'bioflow')
sources_location = os.path.join(storage_location, 'sources')
output_location = os.path.join(storage_location, 'outputs')
internal_storage = os.path.join(storage_location, '.internal')

biogrid_only = False  # if set to true, laplacian and adjacency will use only BioGrid data
