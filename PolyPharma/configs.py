'''
Created on 12 mai 2013

@author: Andrei Kucharavy

Holds all the configurations of the environmental variables for the whole project
'''

# If windows filesystem:
# dbLocation='sqlite:///C:\\Users\\User\\Documents\\UCSD\\DB\\initdb'
# If Skaggs Linux filesystem:
#dbLocation='sqlite:////home/akucahravy/DB/initdb'
# If Linux filesystem on Asus:
dbLocation='sqlite:////media/andrei/OS/Users/Andrei/Documents/UCSD/DB/initdb'

GeneOntology='/home/andrei/workspaces/UCSD/gene_ontology.1_0.obo'
ReactomeBioPax='/home/andrei/workspaces/UCSD/Parsing_Reactome/Homo sapiens.owl'
UNIPROT_text='/home/andrei/workspaces/UCSD/uniprot_sprot.dat'
Hint_csv='/home/andrei/workspaces/UCSD/sapiens_curated-interactome.csv' #attention, for me it is tab-separated

# Targets file is assumed to be tab-separated, with the first column containing fuzzy versions of "names" of "gene names"
# from uniport and the next three - information relative to binding to this target. 
Targets_File='/home/andrei/workspaces/UCSD/NeflanavirSource.csv'

# ExactDict is the dictionnary used to perform a precise matching between the fuzzy target names and the SwissProt IDs required
# for a lookup in the database
from TargetPreProcessing.neflanavir_parser import subdict
Targets_dict = subdict

IDFilter=['5379',     # Small_Molecule_Collection Purine nucleotide
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
               ]
