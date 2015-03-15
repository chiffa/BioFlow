PolyPharma Project
******************

Description:
============

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

The full documentation is available at [readthedocs.org](http://polypharma.readthedocs.org/ RTFD)

Project dependencies:
=====================

**System-level packages:**

* Python 2.7.x (x86_64 )
* Java JDK 1.6/1.7
* Neo4j 1.xx (Note: if you need to install a 2.xx series, please install [gremlin parser engine](https://github.com/neo4j-contrib/gremlin-plugin) too)
* MongoDB
* suitesparese & suitesparse-dev
* LAPACK, ATLAS, BLAS (only if you are installing Scipy manually)


WARNING: all the follwoing is better done within a venv, in order to avoid interfering with other versions of Python
present on the machine


**Packages installed via Pip:**

* Cython
* NumPy (latest x86_64)
* SciPy (latest x86_64)
* bulbs (via pip)
* matplotlib
* python-Levenshtein (via pip)
* pymongo
* requests
* Sphinx (documentation build)
* scikits.sparse (for the cholesky decomposition of a symetric matrix)
* click


Please note that scikits.sparse requires Cython installation via pip and suitesparse-dev(el) installation.

This last step is best done via package manager:

On *Debian*:   ```  $ sudo apt-get install suitesparse suitesparse-dev ```

On *Fedora* / *RHCP* / *CentOS*:    ```  $ sudo yum install suitesparse suitesparse_devel ```


It is recommended to use Gephi for the analysis of the outputs, provided how easy it is for it to 


Required files:
===============
I will write a script to download and install those files automatically from within the application later on,
provided that there are some renaming conventions


**Databases:**
* OBO 1.2 file of GO terms and relations, available at: http://www.geneontology.org/GO.downloads.ontology.shtml
* UNIUPROT-SWISSPROT .txt text database file available at: http://www.uniprot.org/downloads
* Reactome.org "Events in the BioPax level 3" file, available at: http://www.reactome.org/download/index.html
* HiNT binary interaction files for the organisms of interest, availble at: http://hint.yulab.org/batch.html
* BioGRID ALL_ORGANISMS file in the tab2 format, available at http://thebiogrid.org/download.php
* Gene-chromosome mapping files from the Uniprot documentation: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ (needed only for working on aneuploidy)
* Organism-specific protein aboundance files, available at: http://pax-db.org/#!downloads

Please, pay attention to correctly rename the files or edit the files names so that the they correspond to the files specified in sources.ini
file in the `PolyPharma/configs` directory of the application

Basic usage:
============
Start up the neo4j database and the MonogoDB on their default ports. If those ports are not available for some reason, 
please modify the `servers.ini` file in the `/PolyPharma/configs` directory.

**Building the main database**
If you are using the application for the first time on your computer, execute the
```> Python -m PolyPharma.neo4j_Declarations.Import_commander ```
This will build the databases for use with the applications.

** Analysing the mRNA analysis tab file**
Indicate the file to use in the `PolyPharma/configs.py` folder as the RNA_source variable
Configure the expected counts groups and desired intergroup comparisons in the `PolyPharma/PreProcessing/RNA_counts_parser.py` folder 

```> Python -m PolyPharma.PreProcessing.RNA_counts_parser ```

Now, call the auto-analyze routines for the annotation analysis or interactome analysis:
```> Python -m PolyPharma.neo4j_analyzer.knowledge_access_analysis ```

```> Python -m PolyPharma.neo4j_analyzer.interactome_analysis ```

**Command Line interface**
Alternatively, you have direct command line interface available from the root of the project that should avoid you most
of the pain with dealing with inner functions:


*Database creation*

Download the external datasets and build the local master database:

    > python CLUI.py initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver mongodb://localhost:27017/

    > python CLUI.py downloaddbs

    > python CLUI.py setorgconfs --organism [mouse, human, yeast]
    
    > python CLUI.py loadneo4j


*Access to low-level structure

Export the organism-specific interactome or concept-entity relationships as a Python-Scipy sparse matrix object:

    > python CLUI.py extractmatrix --interactome/--annotmap > path to a picke dump of the sparse matrix and name map

Map a list of heterogeneous identifiers to the database-specific ids:

    > python CLUI.py mapids /path/to/my.input.file.tsv > path/to/my.output.file


*Statistical analysis*

 
Analyze a list of genes with an optional background:

    > python CLUI.py analyze --interactome/--annotmap --background /path/to/background.input.file --depth 20 --processors 2 path/to/hits.input.file

The resulting significance data can be seen as the output and the related analyzis .gdf files can be vound in the /outputs folder and 

The full documentation is available at [readthedocs.org](http://polypharma.readthedocs.org/ RTFD)