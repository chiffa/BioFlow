'''
Created on 12 mai 2013
@author: Andrei Kucharavy
Holds all the configurations of the environmental variables for the whole project
'''

from PolyPharma.Utils.ConfigParsers.Configs_parser import parse_configs, sourcefile_compilator
from pprint import PrettyPrinter
from os import path

Servers, Options, Sources, Predictions = parse_configs()

SQLite_location = Servers['PRODUCTION']['local_sqlite']
MongoDB_url = Servers['PRODUCTION']['mongodb_server']
neo4j_server = Servers['PRODUCTION']['server_neo4j']
ReadSourceDBs = sourcefile_compilator(Sources)

GeneOntology = ReadSourceDBs['GO']
ReactomeBioPax = ReadSourceDBs['REACTOME']
UNIPROT_source = ReadSourceDBs['UNIPROT']
Hint_csv = ReadSourceDBs['HINT']  #attention, for me it is tab-separated
Protein_aboundances = ReadSourceDBs['ABOUNDANCES']
sedEffFileName = ReadSourceDBs['SIDER']  #TODO: improve mappings from Drugs to secondary effects.

ReadSourcePredictions = sourcefile_compilator(Predictions)

# Targets file is assumed to be tab-separated, with the first column containing fuzzy versions of "names" of "gene names"
# from uniport and the next three - information relative to binding to this target.
Targets_File = ReadSourcePredictions['NEFLANAVIR']
# File from which to load the names of the 300 most frequent targets
Targets_File2 = ReadSourcePredictions['OVERINGTON']

# ExactDict is the dictionnary used to perform a precise matching between the fuzzy target names and the SwissProt IDs required
# for a lookup in the database
from TargetPreProcessing.neflanavir_parser import subdict
Targets_dict = subdict
from TargetPreProcessing.Overington_parser import subdict2
Targets_dict2 = subdict2

from pymongo import MongoClient
client = MongoClient(MongoDB_url)
db = client.PolyPharma_database
# TODO: see what we are going to do with versionning
ref_coll = db.refrence_v_0_3
data_coll = db.data_v_0_3


# Refers to the groups of links between the nodes that should be treated in the same manner
edge_type_filters = {
    "Group" : ["is_part_of_collection"],                                  # Group relation group
    "Same" : ["is_same"],                                                 # Same relation group
    "Reaction" : ["is_Catalysant", "is_reaction_particpant"],             # Reaction relation group
    "Contact_interaction" : ["is_part_of_complex", "is_Regulant"],        # Contact_interaction relation group
    "HiNT_Contact_interaction" : ["is_interacting"],                      # Contact_interaction relation group
    "possibly_same" : ["is_possibly_same"],
    }


# Coefficients values for the value_Matrix
Adjacency_Martix_Dict = {"Group":0.5,
             "Same":1,
             "Reaction":0.33,
             "Contact_interaction":0.33,
             "possibly_same":0.1,
             }

# Coefficients values for the conductance_Matrix
Conductance_Matrix_Dict = {"Group":0.5,
             "Same":100,
             "Reaction":1,
             "Contact_interaction":1,
             "possibly_same":0.1,
             }

class Dumps(object):
    prefix = path.abspath('./dumps/')
    matrix_LS = 'dump5.dump'
    matrix_corrs = 'dump2.dump'
    eigen_VaMat = 'eigen_valmat.csv'
    eigen_ConMat = 'eigen_conmat.csv'
    val_eigen = 'pickleDump.dump'
    cond_eigen = 'pickleDump_0_5.dump'
    ValMat = 'pickleDump3.dump'
    ConMat = 'pickleDump4.dump'
    UniP_att = 'UP_Attach.dump'
    Main_Connex_group = 'Connex_group.dump'
    Forbidden_IDs = prefix+'ForbiddenIDs.dump'


Leg_ID_Filter = {
    'Reactome:' : ['H+','ATP', 'GTP', 'Pi', 'H2O', 'ADP', 'PPi', 'GDP', 'O2', 'CO2', 'NTP',]
}

# Sh***, those are ID values in a database version, they are not strict



IDFilter = [  '5379',     # Small_Molecule_Collection Purine nucleotide
              '5298',     # Small_Molecule_Collection ADP, GDP, CDP, UDP
              '816',      # Small_Molecule ATP
              '2419',     # Small_Molecule ADP
              '5391',     # Small_Molecule_Collection (d)ADP
              '40020',    # Small_Molecule H2O
              '2047',     # Small_Molecule H+
              '5296',     # Small_Molecule_Collection ADP, GDP, CDP, UDP
              '2795',     # ('Small_Molecule', u'H2O', 'cytosol')  (All the follwong retrieved by the eignevalue approach)
              '2217',     # ('Small_Molecule', u'H2O', 'lysosomal lumen')
              '806',      # ('Small_Molecule', u'Pi', 'cytosol')
              '4055',     # ('Small_Molecule', u'H2O', 'extracellular region')
              '2880',     # ('Small_Molecule', u'CO2', 'cytosol')
              '2317',     # ('Small_Molecule', u'O2', 'cytosol')
              '2803',     # ('Small_Molecule', u'PPi', 'cytosol')
              '2045',     # ('Small_Molecule', u'H+', 'extracellular region')
              '812',      # ('Small_Molecule', u'ATP', 'nucleoplasm')
              '1758',     # ('Small_Molecule', u'ATP', 'mitochondrial matrix')
              '1740',     # ('Small_Molecule', u'H+', 'mitochondrial matrix')
              '992',      # ('Small_Molecule', u'H2O', 'endoplasmic reticulum lumen')
              '2429',     # ('Small_Molecule', u'H2O', 'nucleoplasm')
              '1240',     # ('Small_Molecule', u'O2', 'peroxisomal matrix')
              '818',      # ('Small_Molecule', u'Pi', 'nucleoplasm')
              '810',      # ('Small_Molecule', u'ADP', 'nucleoplasm')
              '2727',     # ('Small_Molecule', u'H2O', 'mitochondrial matrix')
              '1601',     # ('Small_Molecule', u'H2O', 'peroxisomal matrix')
              '2157',     # ('Small_Molecule', u'O2', 'endoplasmic reticulum lumen')
              '2741',     # ('Small_Molecule', u'ADP', 'mitochondrial matrix')
              '2347',     # ('Small_Molecule', u'Pi', 'mitochondrial matrix')
              '1050',     # ('Small_Molecule', u'H+', 'endoplasmic reticulum lumen')
              '2151',     # ('Small_Molecule', u'CO2', 'endoplasmic reticulum lumen')
              '872',      # H+
               ]



if __name__ == "__main__":
    pp=PrettyPrinter(indent=4)
    pp.pprint((Servers,Options,Sources))
