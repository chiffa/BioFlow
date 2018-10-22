|License Type| |Branch status| |Python version|

BioFlow Project
===============

Information Flow Analysis in biological networks

|Build Status| |Coverage Status|  |Duplicate Lines| |Code Health|

Description:
------------

This project's goal is to predict a systemic effect of massive gene
perturbation, whether triggered by a drug, causative mutation or a disease
(such as cancer or disease with complex genetic background). It's main intended
uses are the reduction of high-throughput experiments hit lists, in-silico prediction
of de-novo drug toxicity based on their protein binding profile and retrieval of
most likely pathways explaining a phenotype of interest from a complex genotype.

Its main advantage is the integration of quantitative computational
predictions with prior biological knowledge and ability to integrate
such diverse source of knowledge as databases, simulation, publication
data and expert knowledge.

Unlike similar solutions, it provides several levels of access to the underlying data (integrated
database instance with graph visualization, courtesy of `neo4j graph platform <https://neo4j.com/>`__,
as well as python `numpy <http://www.numpy.org/>`__/`scikits <https://www.scipy.org/>`__
sparse adjacency and laplacian graphs.

The application is currently under development (alpha), hence the API is unstable and can be changed
at any point without notice. If you are using it, please pin the version/commit number. If you
run into issues, please fill the github ticket.

The license is BSD 3-clause, in case of academic usage, please cite the *url* of this repository
(publication is in preparation). The full API documentation is available at
`readthedocs.org <http://bioflow.readthedocs.org/en/latest/>`__.

Installation walk-through:
--------------------------

Ubuntu desktop:
```````````````

1) Install the Anaconda python 2.7 and make it your default python. The full process is explained
`here <https://docs.anaconda.com/anaconda/install/linux/>`__

2) Isnstall libsuitesparse: ::

    > apt-get -y install libsuitesparse-dev

3) Install neo4j: ::

    > wget -O - https://debian.neo4j.org/neotechnology.gpg.key | sudo apt-key add -
    > echo 'deb https://debian.neo4j.org/repo stable/' | sudo tee /etc/apt/sources.list.d/neo4j.list
    > sudo apt-get update
    > sudo apt-get install neo4j

4) Install MongDB (Assuming Linux 18.04 - if not, see
`here <https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/>`__): ::

    > sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 9DA31620334BD75D9DCB49F368818C72E52529D4
    > echo "deb [ arch=amd64 ] https://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/4.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.0.list
    > sudo apt-get update
    > sudo apt-get install -y mongodb-org

For more information, refer to the `installation guide
<http://bioflow.readthedocs.org/en/latest/guide.html#installation-and-requirements>`__

5) Finally, install BioFlow: ::

    > pip install BioFlow


Docker:
```````


Usage walk-through:
-------------------

Fire up the databases ::

    > /home/ank/neo4j-yeast/bin/neo4j start
    > nohup /home/ank/mongodb/bin/mongod &

Setup environment ::

    > bioflow initialize --/home/ank/data_store
    > bioflow downloaddbs
    > bioflow setorg yeast
    > bioflow loadneo4j

For more information about data and config files, refer to the `data and database guide
<http://bioflow.readthedocs.org/en/latest/guide.html#data-and-databases-setup>`__

Set the set of perturbed proteins on which we would want to base our analysis ::

    > bioflow setsource /home/ank/source_data/perturbed_proteins_ids.csv

Build network interfaces ::

    > bioflow extractmatrix --interactome
    > bioflow extractmatrix --annotome

Perform the analysis::

    > bioflow analyze --matrix interactome --depth 24 --processors 4
    > bioflow analyze --matrix annotome --depth 24 --processors 4

The results of analysis will be available in the output folder, and printed out to the standard
output.

For more details or usage as a library, refer to the `usage guide
<http://bioflow.readthedocs.org/en/latest/guide.html#basic-usage>`__.

.. |License Type| image:: https://img.shields.io/badge/license-BSD3-blue.svg
   :target: https://github.com/chiffa/BioFlow/blob/master/License-new_BSD.txt
.. |Build Status| image:: https://travis-ci.org/chiffa/BioFlow.svg?branch=master
   :target: https://travis-ci.org/chiffa/BioFlow
.. |Coverage Status| image:: https://coveralls.io/repos/chiffa/BioFlow/badge.svg?branch=master&service=github
   :target: https://coveralls.io/github/chiffa/BioFlow?branch=master

.. |Duplicate Lines| image:: https://img.shields.io/badge/duplicate%20lines-11.45%25-yellowgreen.svg
   :target: http://clonedigger.sourceforge.net/
.. |Code Health| image:: https://landscape.io/github/chiffa/BioFlow/master/landscape.svg?style=flat
   :target: https://landscape.io/github/chiffa/BioFlow/master

.. |Python version| image:: https://img.shields.io/badge/python-2.7-blue.svg
   :target: https://www.python.org/downloads/release/python-2715/
.. |Branch Status| image:: https://img.shields.io/badge/status-alpha-red.svg
   :target: https://www.python.org/downloads/release/python-2715/
