[![License
Type](https://img.shields.io/badge/license-BSD3-blue.svg)](https://github.com/chiffa/BioFlow/blob/master/License-new_BSD.txt)
[![Python
version](https://img.shields.io/badge/python-2.7-blue.svg)](https://www.python.org/downloads/release/python-2715/)
[![Documentation Status](https://readthedocs.org/projects/bioflow/badge/?version=latest)](https://bioflow.readthedocs.io/en/latest/?badge=latest)
[![Branch
Status](https://img.shields.io/badge/status-alpha-red.svg)](https://www.python.org/downloads/release/python-2715/)

BioFlow Project
===============

Information Flow Analysis in biological networks

[![Build
Status](https://travis-ci.org/chiffa/BioFlow.svg?branch=master)](https://travis-ci.org/chiffa/BioFlow)
[![Coverage
Status](https://coveralls.io/repos/chiffa/BioFlow/badge.svg?branch=master&service=github)](https://coveralls.io/github/chiffa/BioFlow?branch=master)
[![Duplicate
Lines](https://img.shields.io/badge/duplicate%20lines-11.45%25-yellowgreen.svg)](http://clonedigger.sourceforge.net/)
[![Code
Health](https://landscape.io/github/chiffa/BioFlow/master/landscape.svg?style=flat)](https://landscape.io/github/chiffa/BioFlow/master)

Description:
------------

This project's goal is to predict a systemic effect of massive gene
perturbation, whether triggered by a drug, causative mutation or a
disease (such as cancer or disease with complex genetic background).
It's main intended uses are the reduction of high-throughput experiments
hit lists, in-silico prediction of de-novo drug toxicity based on their
protein binding profile and retrieval of most likely pathways explaining
a phenotype of interest from a complex genotype.

Its main advantage is the integration of quantitative computational
predictions with prior biological knowledge and ability to integrate
such diverse source of knowledge as databases, simulation, publication
data and expert knowledge.

Unlike similar solutions, it provides several levels of access to the
underlying data (integrated database instance with graph visualization,
courtesy of [neo4j graph platform](https://neo4j.com/), as well as
python [numpy](http://www.numpy.org/)/[scikits](https://www.scipy.org/)
sparse adjacency and laplacian graphs.

The application is currently under development (alpha), hence the API is
unstable and can be changed at any point without notice. If you are
using it, please pin the version/commit number. If you run into issues,
please fill the github ticket.

The license is BSD 3-clause, in case of academic usage, please cite the
*url* of this repository (publication is in preparation). The full API
documentation is available at
[readthedocs.org](http://bioflow.readthedocs.org/en/latest/).

Installation walk-through:
--------------------------

### Ubuntu desktop:

1)  Install the Anaconda python 2.7 and make it your default python. The
    full process is explained
    [here](https://docs.anaconda.com/anaconda/install/linux/)
2)  Isnstall libsuitesparse:

        > apt-get -y install libsuitesparse-dev

3)  Install neo4j:

        > wget -O - https://debian.neo4j.org/neotechnology.gpg.key | sudo apt-key add -
        > echo 'deb https://debian.neo4j.org/repo stable/' | sudo tee /etc/apt/sources.list.d/neo4j.list
        > sudo apt-get update
        > sudo apt-get install neo4j

4)  Install MongDB (Assuming Linux 18.04 - if not, see
    [here](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/)):

        > sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 9DA31620334BD75D9DCB49F368818C72E52529D4
        > echo "deb [ arch=amd64 ] https://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/4.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.0.list
        > sudo apt-get update
        > sudo apt-get install -y mongodb-org

For more information, refer to the [installation
guide](http://bioflow.readthedocs.org/en/latest/guide.html#installation-and-requirements)

5)  Finally, install BioFlow: :

        > pip install BioFlow

### Docker:

If you want to build locally (notice you need to issue docker commands
with the actual docker-enabled user; usually prepending sudo to the
commands):

    > cd <BioFlow installation folder>
    > docker build -t
    > docker run bioflow
    > docker-compose build
    > docker-compose up -d

If you want to pull from dockerhub or don't have access to BioFlow
installation directory:

    > wget https://github.com/chiffa/BioFlow/blob/master/docker-compose.yml
    > docker-compose build
    > docker-compose up -d

Usage walk-through:
-------------------

### Python scripts:

This is the recommended method for using BioFlow.

Import the minimal dependencies:

    > from bioflow.annotation_network.knowledge_access_analysis import auto_analyze as knowledge_analysis
    > from bioflow.molecular_network.interactome_analysis import auto_analyze as interactome_analysis
    > from bioflow.utils.io_routines import get_source_bulbs_ids
    > from bioflow.utils.top_level import map_and_save_gene_ids, rebuild_the_laplacians

Set static folders and urls for the databases & pull the online
databases:

    > set_folders('~/support') # script restart here is required to properly update all the folders
    > pull_online_dbs()

Set the organism (human, S. Cerevisiae):

    > build_source_config('human')  # script restart here is required to properly update all the folders

Map the hits and the background genes (available through an experimental
technique) to internal IDs:

    > map_and_save_gene_ids('path_to_hits.csv', 'path_to_background.csv')

BioFlow expects the csv files to contain one gene per line. It is
capable of mapping genes based on the following ID types:

-   Gene names
-   HGCN symbols
-   PDB Ids
-   ENSEMBL Ids
-   RefSeq IDs
-   Uniprot IDs
-   Uniprot accession numbers

Preferred mapping approach is through HGCN symbols and Gene names.

Rebuild the laplacians (not required unless background Ids List has been
changed):

    > rebuild_the_laplacians(all_detectable_genes=background_bulbs_ids)

Launch the analysis itself for the information flow in the interactome:

    > interactome_analysis([hits_ids],
                          desired_depth=9,
                          processors=3,
                          background_list=background_bulbs_ids,
                          skip_sampling=False,
                          from_memoization=False)

Launch the analysis itself for the information flow in the annotation
network (experimental):

    > knowledge_analysis([hits_ids],
                        desired_depth=20,
                        processors=3,
                        skip_sampling=False)

Where:

hits\_ids
:   list of hits

desired\_depth
:   how many samples we would like to generate to compare against

processors
:   how many threads we would like to launch in parallel (in general 3/4
    works best)

background\_list
:   list of background Ids

skip\_sampling
:   if true, skips the sampling of background set and retrieves stored
    ones instead

from\_memoization
:   if true, assumes the information flow for the hits sample has
    already been computed

BioFlow will print progress to the StdErr from then on and will output
to the user's home directory, in a folder called 'outputs\_YYYY-MM\_DD
\<launch time\>':

-   .gdf file with the flow network and relevance statistics
    (Interactome\_Analysis\_output.gdf)
-   visualisation of information flow through nodes in the null vs hits
    sets based on the node degree
-   list of strongest hits (interactome\_stats.tsv)

The .gdf file can be further analysed with more appropriate tools.

### Command line:

Setup environment (likely to take a while top pull all the online
databases): :

    > bioflow initialize --~/data_store
    > bioflow downloaddbs
    > bioflow setorg human
    > bioflow loadneo4j

For more information about data and config files, refer to the [data and
database
guide](http://bioflow.readthedocs.org/en/latest/guide.html#data-and-databases-setup)

Set the set of perturbed proteins on which we would want to base our
analysis :

    > bioflow setsource /home/ank/source_data/perturbed_proteins_ids.csv

Build network interfaces :

    > bioflow extractmatrix --interactome
    > bioflow extractmatrix --annotome

Perform the analysis:

    > bioflow analyze --matrix interactome --depth 24 --processors 4
    > bioflow analyze --matrix annotome --depth 24 --processors 4

The results of analysis will be available in the output folder, and
printed out to the standard output.

### Post-processing:

The .gdf file format is one of the standard format for graph exchange.
It contains the following columns for the nodes:

-   node ID
-   information current passing through the node
-   node type
-   legacy\_id (most likely Uniprot ID)
-   degree of the node
-   whether it is present or not in the hits list (source)
-   p-value, comparing the information flow through the node to the flow
    expected for the random set of genes
-   -log10(p\_value) (p\_p-value)
-   rel\_value (information flow relative to the flow expected for a
    random set of genes)
-   std\_diff (how many standard deviations above the flow for a random
    set of genes the flow from a hits list is)

The most common pipleine involves using [Gephi open graph visualization
platform](https://gephi.org/):

-   Load the gdf file into gephy
-   Filter out all the nodes with information flow below 0.05 (Filters
    \> Atrributes \> Range \> current)
-   Perform clustering (Statistics \> Modularity \> Randomize & use
    weights)
-   Filter out all the nodes below a significance threshold (Filters \>
    Attributes \> Range \> p-value)
-   Set Color nodes based on the Modularity Class (Nodes \> Colors \>
    Partition \> Modularity Class)
-   Set node size based on p\_p-value (Nodes \> Size \> Ranking \>
    p\_p-value )
-   Set text color based on whether the node is in the hits list (Nodes
    \> Text Color \> Partition \> source)
-   Set text size based on p\_p-value (Nodes \> Text Size \> Ranking \>
    p\_p-value)
-   Show the lables (T on the bottom left)
-   Set labes to the legacy IDs (Notepad on the bottom)
-   Perform a ForeAtlas Node Separation (Layout \> Force Atlas 2 \>
    Dissuade Hubs & Prevent Overlap)
-   Adjust label size
-   Adjust labels position (Layout \> LabelAdjust)

For more details or usage as a library, refer to the [usage
guide](http://bioflow.readthedocs.org/en/latest/guide.html#basic-usage).
