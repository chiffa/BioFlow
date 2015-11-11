[![License Type](https://img.shields.io/badge/license-BSD3-blue.svg)](https://github.com/chiffa/Karyotype_retriever/blob/master/License-BSD3)
[![Build Status](https://travis-ci.org/chiffa/BioFlow.svg?branch=v0.03)](https://travis-ci.org/chiffa/BioFlow)
[![Coverage Status](https://coveralls.io/repos/chiffa/PolyPharma/badge.svg?branch=master&service=github)](https://coveralls.io/github/chiffa/PolyPharma?branch=master)
# BioFlow Project
PInformation Flow Analysis toolkit for biological networks

## Description:

This code is currently under development and thus it isn't always stable, thoroughly tested or well documented.

This project's goal is to predict systemic effect of multiple gene perturbation, wherther trigered by a drug or by
a disease (such as cancer or disease with complex genetic background). It's main intended uses are prediction of
drug toxicity of de-novo drugs due to a distributed off-target effect and linkage between a phenotype and a complex
genotype.

It's main advantage is integration of quantitative computational predictions with prior biological knowledge and
ability to integrate such diverse source of knowledge as databases, simulation, publication data (currently in dev)
and expert knowledge.

This application is provided as-is, without any warranties or support. Use it at your onw risk.

However, if you are willing to test it and encounter problems or are willing to provide feedback, please fill in
an issue ticket on GitHub and I will be glad to assist you in the measure of my possibilities.

The license is BSD, but in case of academic usage, please cite the *url* (publication is in preparation).

The full documentation is available at [readthedocs.org](http://bioflow.readthedocs.org/ RTFD)

## Installation:

Before you start using the project, you will need to install several dependencies on which the project relies in order 
to function properly.

Currently, only deployment on Linux Debian/RHCP/CentOS/Fedora is supported because of dependencies on scientific
computation packages

**System-level packages:**

* Python 2.7.x (x86_64 )
* Java JDK 1.6/1.7
* Neo4j 1.xx You can download it from [the official Neo4j website](http://neo4j.com/download/other-releases/).
Note: if you want to install a 2.xx series you will also need a [gremlin parser engine](https://github.com/neo4j-contrib/gremlin-plugin).
In case your neo4j version fails to with a 'java heap exception', please considers installing a [frozen 1.9.6 community neo4j version](https://github.com/chiffa/neo4j-community-1.9.6)
* MongoDB ([installation instructions here](https://docs.mongodb.org/manual/administration/install-on-linux/))
* LAPACK, ATLAS, BLAS (only if you are installing Scipy manually)
* suitesparese & suitesparse-dev (needed for scikits.sparse compilation. if you are installing a wheel of scikits.sparse, you don't need them)

This last step is best done via package manager:

On *Debian*:   ```  $ sudo apt-get install suitesparse suitesparse-dev ```

On *Fedora* / *RHCP* / *CentOS*:    ```  $ sudo yum install suitesparse suitesparse_devel ```


**Packages installed via Pip:**

All the Python packages are better off stored in a, so they don't interfere with other versions of Python present on the
 machine. For a more simple installation, you can install [Anaconda Python Distribution](https://www.continuum.io/downloads),
 which will install the following packages:

* Cython
* NumPy (latest x86_64)
* SciPy (latest x86_64)
* matplotlib
* Sphinx (documentation build)

After which, you can `pip install` the remaining packages

* bulbs 
* python-Levenshtein
* pymongo
* requests
* scikits.sparse
* click


Finally, I would recommend using [Gephi](http://gephi.github.io/users/download/) to analyse the output .gdf graphs.
 It is very easy to use, well-supported and possess a large library of plug-ins.


## Datasets used as backbone:

* OBO 1.2 file of GO terms and relations, available at: http://www.geneontology.org/GO.downloads.ontology.shtml
* UNIUPROT-SWISSPROT .txt text database file available at: http://www.uniprot.org/downloads
* Reactome.org "Events in the BioPax level 3" file, available at: http://www.reactome.org/download/index.html
* HiNT binary interaction files for the organisms of interest, availble at: http://hint.yulab.org/batch.html
* BioGRID ALL_ORGANISMS file in the tab2 format, available at http://thebiogrid.org/download.php
* Gene-chromosome mapping files from the Uniprot documentation: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ (needed only for working on aneuploidy)
* Organism-specific protein aboundance files, available at: http://pax-db.org/#!downloads

You have a choice of either donwloading and these databases and installing them manually with respect to 

## Basic usage:

Start up the neo4j database and the MonogoDB on their default ports. If those ports are not available for some reason, 
please modify the `servers.ini` file in the `/PolyPharma/configs` directory. If you are loading a particularly large
dataset into neo4j, you might need to adjust the configurations with which neo4j is launched

You have a choice of using either Python binding of methods or a command line interface. If you installed this package with
pip, your will be able to call the command line interfaces with

    > truegrid command args
    
Otherwise, you can access them with
    
    > python CLUI.py command args
    
when you are executing the command lines by using Python bindings, don't forget to add the `-m` argument between
`Python` and the location of the command you are issuing.

    > Python -m PolyPharma.module.command


### Downloading the datasets:

*Database creation*

Create the configuration files containing the data:

    > python CLUI.py initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver mongodb://localhost:27017/
    
    > python -m PolyPharma.Utils.ConfigsIO.set_folders()
    
Download the databases:

    > python CLUI.py downloaddbs
    
    > python -m PolyPharma.Utils.ConfigsIO.StructureGenerator.pull_online_DBs()

For now, the syustem will download all the required files, then fail when trying to download 'ABOUNDANCE' file class. 
    
Create the proper configuration file for the desired organism

    > python CLUI.py setorgconfs --organism [mouse, human, yeast]
    
    > python - m PolyPharma.Utils.ConfigsIO.build_source_config('yeast')
    
Alternatively all of the above can be executed (for yeast),

    > python- m PolyPharma.Utils.ConfigsIO 
    
Provided that Uniprot.dat is a rather big file (~3 Gb as of late 2015), it might get broken on the download and you
might want to check that it is a correct size


### Building the database:

If you are using the application for the first time on your computer, you will need to load all the data that is
contained in the datastore files you've donwloaded previously and cross-reference them

    > Python -m PolyPharma.neo4j_Importers.Import_commander

    > python CLUI.py loadneo4j


### Accessing low-level structure of the interactome:

Export the organism-specific interactome or concept-entity relationships as a Python-Scipy sparse matrix object:

    > python CLUI.py extractmatrix --interactome/--annotmap > path to a picke dump of the sparse matrix and name map

    > python -m 

Map a list of heterogeneous identifiers to the database-specific ids:

    > python CLUI.py mapids /path/to/my.input.file.tsv > path/to/my.output.file
    
    > python -m
    

### High-level analysis:

Indicate the file to use in the `PolyPharma/configs.py` folder as the RNA_source variable
Configure the expected counts groups and desired intergroup comparisons in the `PolyPharma/PreProcessing/RNA_counts_parser.py` folder 

    > Python -m PolyPharma.PreProcessing.RNA_counts_parser

Now, call the auto-analyze routines for the annotation analysis or interactome analysis:

    > Python -m PolyPharma.neo4j_analyzer.knowledge_access_analysis

    > Python -m PolyPharma.neo4j_analyzer.interactome_analysis

 
Analyze a list of genes with an optional background:

    > python CLUI.py analyze --interactome/--annotmap --background /path/to/background.input.file --depth 20 --processors 2 path/to/hits.input.file

The resulting significance data can be seen as the output and the related analyzis .gdf files can be found in the /outputs folder.


Full API documentation of underlying libraries is available at [readthedocs.org](http://polypharma.readthedocs.org/ RTFD)


## Future developments:

 - Automatic change of neo4j instances for different organisms (currently needs to be done manually)
 