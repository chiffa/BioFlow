Installation and requirements:
==============================

Right now the installation is possible only on the Debian and RHCP flavors of Linux because of
dependencies on the scientific computation libraries. If you desire to install on other OSes, as
long as you are able to satisfy the Neo4j/Mongo requirements and properly install all the Python
modules, the library should work fine.

If you use Ubuntu 14.04 or are running a docker image of Ubuntu 14.04, you can just execute the
``ubuntu_14_04_setup.sh`` script: ::

    > sh docker_set_up.sh

Otherwise, you will need to perform two-folded set-up procedure:

Install the system-level packages:
----------------------------------

- Python 2.7
- Oracle Java JDK 1.7
- Neo4j 1.9.6 (latest bugless release containing gremlin). You can get a frozen version `here
<https://github.com/chiffa/neo4j-community-1.9.6>`__
- Mongo database (installation instructions
   `here <https://docs.mongodb.org/manual/administration/install-on-linux/>`__)
- LAPACK, ATLAS and BLAS (if you will be installing Numpy and Scipy manually)
- (lib)suitesparse-dev (needed for scikits.sparse
   compilation)

If you are operating on debian (here ubuntu 14.04), you can refer to
``docker_set_up.sh`` for the instructions that are issued to get each part

Install python libraries:
-------------------------

If you install this package via ``pip install ``, all the dependencies should be installed.
Alternatively, you can install them by issuing ::

    > pip install requirements -r requirements.txt

If you wish to install all the dependencies manually, you will need:

    -  Cython
    -  NumPy (latest x86\_64)
    -  SciPy (latest x86\_64)
    -  matplotlib
    -  Sphinx (documentation build)
    -  bulbs
    -  python-Levenshtein
    -  pymongo
    -  requests
    -  scikits.sparse
    -  click

Finally, if you want to install them quickly or just isolate from the rest of your Python
environment, `conda <https://www.continuum.io/downloads>`__ virtual environment provide a very good
 way of doing it:

Detailed instructions for ubuntu 14.04:
---------------------------------------

Install basic dependencies for Python modules::

    apt-get -y install python
    apt-get -y install lsof libsm6 libxrender1 libfontconfig1 libglib2.0-0
    apt-get -y install build-essential
    apt-get -y install libsuitesparse-dev

    # in case you are on a docker image of ubuntu 14.04:
    apt-get -y install wget git curl nohup

Install and activate Oracle Java 1.7::

    apt-get -y install software-properties-common
    apt-get -y install python-software-properties
    add-apt-repository ppa:webupd8team/java
    apt-get -y update

    apt-get -y install oracle-java7-installer
    apt-get -y install oracle-java7-set-default

Install neo4j::

    git clone https://github.com/chiffa/neo4j-community-1.9.6.git
    mv neo4j-community-1.9.6 neo4j-yeast

Install mongodb::

    curl -O https://fastdl.mongodb.org/linux/mongodb-linux-x86_64-3.0.8.tgz
    tar -zxvf mongodb-linux-x86_64-3.0.8.tgz
    rm mongodb-linux-x86_64-3.0.8.tgz
    mv mongodb-linux-x86_64-3.0.8 mongodb
    mkdir -p /data/db

Create and activate conda test environments::

    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p /home/ank/miniconda
    export PATH="/home/ank/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    rm miniconda.sh

    conda create -q -n test-environement python="2.7" numpy scipy matplotlib
    source activate test-environement
    conda install python="2.7" cython scikit-learn




Software for graph visualization:
---------------------------------

Network analysis results will be output as `.gdf` files. In order to visualize them, I would
recommend using `Gephi <https://gephi.org/users/download/>`__.


Data and databases setup:
=========================

Assembling the files required for the database creation:
--------------------------------------------------------

In order to build the database, the program is going to look for the following files specified
in the following locations within the PolyPharma/configs/sources.ini::

    * OBO 1.2 file of GO terms and relations, available at: http://www.geneontology.org/GO.downloads.ontology.html
    * will look for at the path in [GO] - "location"

    * UNIUPROT-SWISSPROT .txt text database file available at: http://www.uniprot.org/downloads
    * will look for the files at the path [UNIPROT] - "location"
    * will load the information for the organism with specified NCBI taxonomy identifier from "tax_id"

    * Reactome.org "Events in the BioPax level 3" file, available at: http://www.reactome.org/download/index.html
    * will look for the files at [REACTOME] - "location"
    * will only load the file specified by the "load" parameter

    * HiNT binary interaction files for the organisms of interest, availble at: http://hint.yulab.org/batch.html
    * will look for the files at the path [HINT] - "location"
    * will load the information for the organism with specified NCBI taxonomy identifier from "load" parameter

    * BioGRID ALL_ORGANISMS file in the tab2 format, available at http://thebiogrid.org/download.php
    * will look for for the files at the path [BIOGIRD] - "location"
    * will load only the file specified in the "load" parameter

    * Gene-chromosome mapping files from the Uniprot documentation: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ (needed only for working on aneuploidy)
    * Configured in the [CHROMOSOMES] section
    * currently not ready for public use

    * Organism-specific protein aboundance files, available at: http://pax-db.org/#!downloads
    * Configured in the [ABOUNDACES] section
    * currently not ready for public use


It is possible to specify the file locations and identifiers manually, and then download and install them
to the specified locations manually.

However the following command should be able to do it for you for three commonly used organism (human, mouse, saccharomyces cerevisae),
provided you follow the instructions properly::

    > python CLUI.py initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver mongodb://localhost:27017/

    > python CLUI.py downloaddbs

    > python CLUI.py setorgconfs --organism [mouse, human, yeast]

Typical sources.ini configfile:
-------------------------------

Here is what a typical configfile would look like::

    [REACTOME]
    location = myfolder/External_DBs_Store/Reactome
    load = Mus musculus.owl

    [UNIPROT]
    location = myfolder/External_DBs_Store/Uniprot/uniprot_sprot.dat
    tax_ids = 10090,

    [HINT]
    location = myfolder/E/External_DBs_Store/HiNT
    load = MouseBinaryHQ.txt

    [GO]
    location = myfolder/E/External_DBs_Store/GO/go.obo

    [BIOGIRD]
    location = myfolder/E/External_DBs_Store/BioGIRD
    load = Mus_musculus.tsv

    [CHROMOSOMES]
    location = myfolder/E/External_DBs_Store/Chr_mappings
    load = mouse
    namepattern = mouse

    [ABOUNDANCES]
    location = /home/ank/Documents/External_DBs_Store/Protein_aboundances
    load = 10090

The data relative to the following parameters::

    [REACTOME]

    [UNIPROT]

    [HINT]

    [GO]

    [BIOGIRD]

is critical for any application and must be properly configured and is critical for any application
of the method.

On the other hand the following parameters are here for legacy application reasons and are not currently
documented::

    [CHROMOSOMES]

    [ABOUNDANCES]

the "load" parameter in the "[UNIPROT]" folder requires a trailing comma and can take in multiple arguments
separaged by a comma and a space, in case UNIPROT identifiers of proteins from several organisms are desired
(for instance when host-disease proteome interactions are investigated)

The configuration files might be declared and switchedmanually (only the "source.ini" one will be parsed,
folders such as "sources_organism.ini" will be ignored and can be renamed to "source.ini" quite easily)

It is possible for the users to generate source.ini file for three organisms with the following command::

    python CLUI.py setorgconfs --organism [mouse, human, yeast]

This allows to switch rapidly between different investigated organism.

Please don't forget to switch or purge neo4j databases between organisms, because each organism needs it's own neo4j instance.

Basic usage:
============

Neo4j and mongodb startup:
--------------------------

Start up the neo4j database and the MonogoDB on their default ports. If
those ports are not available for some reason, please modify the
``servers.ini`` file in the ``/PolyPharma/configs`` directory. If you
are loading a particularly large dataset into neo4j, you will need to adjust the java heap space
used to launch the jvm. You can do it by editing the ``$NEO4J_HOME/conf/neo4j-wrapper.conf``
file, and uncommenting + editing the ``wrapper.java.initmemory`` and ``wrapper.java.maxmemory``

Neo4j out of memory error:
--------------------------

In case you are going to work with organisms with large proteomes (mouse, human), neo4j might run
 out of memory and prompt to be restarted with a larger allocation of RAM. In order to correct
 this error, please umncomment and modify the following lines in the ``neo4j-wrapper.conf`` file in
 your neo4j installation instance to  increase the initial and maximum memory for java process
 running neo4j: ::

    wrapper.java.initmemory=16
    wrapper.java.maxmemory=64

Command line:
-------------

An example of usage of the command line interface is given in the Readme, however we will
implement it again here:

Provide local datastore location ::

    > bioflow initialize --/home/ank/data_store

Donwload data repositories to the local datastore ::

    > bioflow downloaddbs

Set organism we want to analyse to yeast ::

    > bioflow setorg yeast

Load the data from the local datastore into the neo4j instance ::

    > bioflow loadneo4j

Set the set of perturbed proteins on which we would want to base our analysis ::

    > bioflow setsource /home/ank/source_data/perturbed_proteins_ids.csv

Buid interactome interface ::

    > bioflow extractmatrix --interactome

Build annotome interface ::

    > bioflow extractmatrix --annotome

Peform the analysis of the set of interest against the matched random sample of size 24, sampled
on 4 processors with respect to the interactome structure ::

    > bioflow analyze --matrix interactome --depth 24 --processors 4

Perform the analysis of the set of interest against the matched random sample of size 24, sampled
 on 4 processors with respect to the annotation structure ::

    > bioflow analyze --matrix annotome --depth 24 --processors 4

The resulst of analysis will be available in the output folder, and printed out to the standard
output.

As a library:
-------------

An example of usage of the library is given in the file called "analysis_pipeline_example.py". To
 rapidly get over it, here is the minimal analysis pipeline:

Setting static folders and urls for the databases ::

    bioflow.configs_manager.set_folders('/home/ank/data_repository',
                                        'http://localhost:7474',
                                        'mongodb://localhost:27017/')

Pulling the online databases ::

    bioflow.configs_manager.StructureGenerator.pull_online_dbs()

Setting the organism to yeast ::

    bioflow.configs_manager.StructureGenerator.build_source_config('yeast')

Clearing the database, if required ::

    bioflow.db_importers.import_main.destroy_db()

Building the neo4j database for a new organism ::

    bioflow.db_importers.import_main.build_db()

Building the interactome interface object ::

    from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=False)
    local_matrix.full_rebuild()

Setting up the reference parameter set for the analysis of annotome ::

    annotation_type = ['biological_process']
    background_set = local_matrix.all_uniprots_bulbs_id_list
    ref_param_set = [['biological_process'], background_set, (1, 1), True, 3]

Building the annotome interface object ::

    from bioflow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface
    annot_matrix = GeneOntologyInterface(*ref_param_set)
    annot_matrix.rebuild()
    annot_matrix.store()

Set the source file of the ids of perturbed proteins ::

    bioflow.neo4j_db.db_io_routines.cast_analysis_set_to_bulbs_ids(
        "/path/to/perturbed/prots_ids.csv")

Get the bulbs ids oif the nodes we would like to analyze ::

    from bioflow.molecular_network.interactome_analysis \
        import auto_analyze as interactome_analysis, get_source_bulbs_ids
    source_bulbs_ids = get_source_bulbs_ids()

Perform the interactome analysis::

    interactome_analysis([source_bulbs_ids], desired_depth=24, processors=6)

Perform the knowledge analysis ::

    from bioflow.annotation_network.knowledge_access_analysis \
        import auto_analyze as knowledge_analysis
    knowledge_analysis([source_bulbs_ids], desired_depth=24, processors=6)


Graph Analysis with Gephi:
--------------------------

The standard usage pipeline for me is to import the .gdf file into `Gephi <https://gephi.org/users/download/>`__ and:

 - Add a filter for current, set it between 2.5% of max and max

 - In statistics, run modularity analysis

 - Set the color of the nodes as Partition/Modularity Class

 - Set the node size as Ranking/rel_value (current relative to what would have been expected in null conditions)

 - Set the label color as Partition/Source (whether the node was in the hits set or not)

 - Set the label size as Ranking/rel_value

 - Add a Filter for p_value, below 0.05

 - Run Force Atlas, if needed enabling the "no overlap" and "dissuade hubs" options

 - Turn on the labels and switch them to "legacy IDs" or "names"