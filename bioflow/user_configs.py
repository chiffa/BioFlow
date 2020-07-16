import os

storage_location = os.path.join(os.environ['HOME'], 'bioflow')
sources_location = os.path.join(storage_location, 'sources')
output_location = os.path.join(storage_location, 'outputs')
internal_storage = os.path.join(storage_location, '.internal')

molecular_edge_filter = ['All',
                         'Reactome',
                         'BioGRID',
                         'HiNT',
                         'TRRUST']

# TODO: this file needs to be saved in the output folder for every run.