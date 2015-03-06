__author__ = 'ank'

from pycallgraph import PyCallGraph, Config, GlobbingFilter
from pycallgraph.output import GephiOutput
import os
from os.path import join

config = Config()

# config.trace_filter = GlobbingFilter(exclude=[
#     'pycallgraph.',
#     '*.__*',
#     'os.*',
#     'sys.*',
#     'json.*',
#     'atexit.*',
#     'calendar.*',
#     'collections.*',
#     'ctypes.*',
#     'pickle.*',
#     'pprint.*',
#     'random.*',
#     'socket.*',
#     'pymongo.*',
#     'threading.*',
#     'urllib.*',
#
#         ])

config.trace_filter = GlobbingFilter(include=['PolyPharma.*'])

lname = os.path.dirname(os.path.realpath(__file__))

gephi = GephiOutput(output_file=join(lname, 'callgraph.gdf'))


with PyCallGraph(output=gephi, config=config):
    # import PolyPharma.configs as mycode
    # import PolyPharma.PreProcessing.RNA_counts_parser
    import PolyPharma.neo4j_Importers.Import_commander as cmder
    cmder.run_diagnostics({})


