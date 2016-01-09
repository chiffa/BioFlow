
|License Type| |Build Status| |Coverage Status| |Branch status| |Code
Issues| |Duplicate Lines|

BioFlow Project
===============

Information Flow Analysis in biological networks

Branch v.0.03 - Unstable

Description:
------------

This project's goal is to predict systemic effect of multiple gene
perturbation, wherther trigered by a drug or by a disease (such as
cancer or disease with complex genetic background). It's main intended
uses are prediction of drug toxicity of de-novo drugs due to a
distributed off-target effect and linkage between a phenotype and a
complex genotype.

It's main advantage is integration of quantitative computational
predictions with prior biological knowledge and ability to integrate
such diverse source of knowledge as databases, simulation, publication
data (currently in dev) and expert knowledge.

The application is currently under development and in alpha stage. However, if you desire to use
it, you are welcome to do so and fill in the tickets if you encounter any issues

The license is BSD, but in case of academic usage, please cite this *url*
(I am currently writing a publication). The full API documentation is available at
`readthedocs.org <http://bioflow.readthedocs.org/en/latest/>`__.

Installation:
-------------

Before you start using the project, you will need to install several
dependencies on which the project relies in order to function properly.

Quick:
``````
If you are using Ubuntu 14.04, you can just run
::

    > sh docker_set_up.sh

Detailed:
`````````
Only deployment on Linux Debian/RHCP/CentOS/Fedora is
supported because of dependencies on scientific computation packages

You will need the following **system-level packages:**

-  Python 2.7.x (x86\_64 )
-  Java JDK 1.6/1.7
-  Neo4j 1.xx You can download it from `the official Neo4j
   website <http://neo4j.com/download/other-releases/>`__. Note: if you
   want to install a 2.xx series you will also need a `gremlin parser
   engine <https://github.com/neo4j-contrib/gremlin-plugin>`__. In case
   your neo4j version fails to with a 'java heap exception', please
   considers installing a `frozen 1.9.6 community neo4j
   version <https://github.com/chiffa/neo4j-community-1.9.6>`__
-  MongoDB (`installation instructions
   here <https://docs.mongodb.org/manual/administration/install-on-linux/>`__)
-  LAPACK, ATLAS, BLAS (only if you are installing Scipy manually)
-  suitesparese & suitesparse-dev (needed for scikits.sparse
   compilation. if you are installing a wheel of scikits.sparse, you
   don't need them)

This last step is best done via package manager:

On *Debian*:

::

    $ sudo apt-get install suitesparse suitesparse-dev``

On *Fedora* / *RHCP* / *CentOS*:

::

    $ sudo yum install suitesparse suitesparse_devel``


You will also need the following **Python packages:**

All the Python packages are better off stored in a virtual environment, so they don't
interfere with other versions of Python present on the machine. For a
more simple installation, you can install `Anaconda Python
Distribution <https://www.continuum.io/downloads>`__, which will install
the following packages:

-  Cython
-  NumPy (latest x86\_64)
-  SciPy (latest x86\_64)
-  matplotlib
-  Sphinx (documentation build)

After which, you can run ``pip install`` the remaining packages

-  bulbs
-  python-Levenshtein
-  pymongo
-  requests
-  scikits.sparse
-  click

Finally, I would recommend using
`Gephi <http://gephi.github.io/users/download/>`__ to analyse the output
.gdf graphs. It is fairly intuitive and easy to use. I am preparing x-networks integration, but
it is still quite far on the desired features list

Datasets used as backbone:
--------------------------
The method relies on a neo4j database hypergraph where it dumps and cross-references several
biological repositories. normally the files will be downloaded and parsed automatically, but if
for some reason the download breaks, you will need to download them manually and add their paths
in ``configs/sources.ini``

-  `OBO 1.2 file of GO terms and relations
    <http://www.geneontology.org/GO.downloads.ontology.html>`__
-  `UNIUPROT-SWISSPROT swissprot.txt text database file
   <http://www.uniprot.org/downloads>`__
-  `Reactome.org "Events in the BioPax level 3" file
   <http://www.reactome.org/download/index.html>`__
-  `HiNT binary interaction files for the organisms of interest
   <http://hint.yulab.org/batch.html>`__
-  `BioGRID ALL\_ORGANISMS file in the tab2 format
   <http://thebiogrid.org/download.php>`__

Basic usage:
------------

Neo4j and mongodb startup:
``````````````````````````

Start up the neo4j database and the MonogoDB on their default ports. If
those ports are not available for some reason, please modify the
``servers.ini`` file in the ``/PolyPharma/configs`` directory. If you
are loading a particularly large dataset into neo4j, you will need to adjust the java heap space
used to launch the jvm. You can do it by editing the ``$NEO4J_HOME/conf/neo4j-wrapper.conf``
file, and uncommenting + editing the ``wrapper.java.initmemory`` and ``wrapper.java.maxmemory``

Command line usage:
```````````````````
You have a choice of using either Python binding of methods or a command
line interface. If you installed this package with pip, your will be
able to call the command line interfaces with

::

    > bioflow command args

Otherwise, you can access them with

::

    > python CLUI.py command args

Execution of a module:
``````````````````````
when you are executing the command lines by using Python bindings, don't
forget to add the ``-m`` argument between ``Python`` and the location of
the command you are issuing.

::

    > Python -m PolyPharma.module.command

Usage of library:
`````````````````
Most of the modules contain the data

step-by-step command line usage:
````````````````````````````````

Downloading the datasets:
~~~~~~~~~~~~~~~~~~~~~~~~~

*Database creation*

Create the configuration files containing the data:

::

    > bioflow initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver
    mongodb://localhost:27017/

Download the databases:

::

    > biogrid downloaddbs

For now, the syustem will download all the required files.

Create the proper configuration file for the desired organism

::

    > biogrid setorgconfs --organism [mouse, human, yeast]

Alternatively all of the above can be executed (for yeast),

::

    > python- m PolyPharma.Utils.ConfigsIO 

Provided that Uniprot.dat is a rather big file (~3 Gb as of late 2015),
it might get broken on the download and you might want to check that it
is a correct size

Building the database:
~~~~~~~~~~~~~~~~~~~~~~

If you are using the application for the first for an organism,
you will need to load all the data that is contained in the datastore
files you've donwloaded previously and cross-reference them

::

    > biogrid loadneo4j

Accessing low-level structure of the interactome:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Export the organism-specific interactome or concept-entity relationships
as a Python-Scipy sparse matrix object:

::

    > biogrid CLUI.py extractmatrix --interactome/--annotmap > path to a picke dump of the sparse
    matrix and name map

Map a list of heterogeneous identifiers to the database-specific ids:

::

    > biogrid mapids /path/to/my.input.file.tsv > path/to/my.output.file

High-level analysis:
~~~~~~~~~~~~~~~~~~~~

Indicate the file to use in the ``PolyPharma/configs.py`` folder as the
RNA\_source variable Configure the expected counts groups and desired
intergroup comparisons in the
``PolyPharma/PreProcessing/RNA_counts_parser.py`` folder

::

    > Python -m PolyPharma.PreProcessing.RNA_counts_parser

Now, call the auto-analyze routines for the annotation analysis or
interactome analysis:

::

    > Python -m PolyPharma.neo4j_analyzer.knowledge_access_analysis

    > Python -m PolyPharma.neo4j_analyzer.interactome_analysis

Analyze a list of genes with an optional background:

::

    > biogrid analyze --interactome/--annotmap --background /path/to/background.input.file --depth
     20 --processors 2 path/to/hits.input.file

The resulting significance data can be seen as the output and the
related analyzis .gdf files can be found in the /outputs folder. statistic analysis will be
printed to the stdout.

Full API documentation of underlying libraries is available at
`readthedocs.org <http://bioflow.readthedocs.org/en/latest/>`__

Future developments:
--------------------

Please see the developper log below!

.. |License Type| image:: https://img.shields.io/badge/license-BSD3-blue.svg
   :target: https://github.com/chiffa/BioFlow/blob/master/License-new_BSD.txt
.. |Build Status| image:: https://travis-ci.org/chiffa/BioFlow.svg?branch=master
   :target: https://travis-ci.org/chiffa/BioFlow
.. |Coverage Status| image:: https://coveralls.io/repos/chiffa/BioFlow/badge.svg?branch=master&service=github
   :target: https://coveralls.io/github/chiffa/BioFlow?branch=master
.. |Branch status| image:: https://img.shields.io/badge/branch_status-refactoring-red.svg
   :target: https://github.com/chiffa/BioFlow/blob/master/README.md
.. |Code Issues| image:: https://www.quantifiedcode.com/api/v1/project/1c3f8cd001a44319abddab249101b646/badge.svg
   :target: https://www.quantifiedcode.com/app/project/1c3f8cd001a44319abddab249101b646
.. |Duplicate Lines| image:: https://img.shields.io/badge/duplicate%20lines-17.66%25-yellow.svg
   :target: http://clonedigger.sourceforge.net/