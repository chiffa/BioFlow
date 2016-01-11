
|License Type| |Build Status| |Coverage Status| |Branch status| |Code
Issues| |Duplicate Lines|

BioFlow Project
===============

Information Flow Analysis in biological networks

Branch v.0.03 - Release candidate for v0.4.0

Description:
------------

This project's goal is to predict systemic effect of multiple gene
perturbation, whether triggered by a drug or by a disease (such as
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

The license is BSD, but in case of academic usage, please cite the *url* of this repository
(I am currently writing a publication). The full API documentation is available at
`readthedocs.org <http://bioflow.readthedocs.org/en/latest/>`__.

Installation walk-through:
--------------------------

If you are on Ubuntu 14.04: ::

    > sh ubuntu_14_04_setup.sh
    > pip install git+https://github.com/chiffa/BioFlow.git

For more information, refer to the `installation guide <http://bioflow.readthedocs
.org/en/latest/guide.html#installation-and-requirements>`__

Usage walk-through:
-------------------

Fire up the databases ::

    > /home/ank/neo4j-yeast/bin/neo4j start
    > nohup /home/ank/mongodb/bin/mongod &

Setup environement ::

    > bioflow initialize --/home/ank/data_store
    > bioflow downloaddbs
    > bioflow setorg yeast
    > bioflow loadneo4j

For more information about data and config files, refer to the `data and database guide
<http://bioflow.readthedocs.org/en/latest/guide.html#data-and-database-setup>`__

Set the set of perturbed proteins on which we would want to base our analysis ::

    > bioflow setsource /home/ank/source_data/perturbed_proteins_ids.csv

Build netowrk interfaces ::

    > bioflow extractmatrix --interactome
    > bioflow extractmatrix --annotome

Peform the analysis::

    > bioflow analyze --matrix interactome --depth 24 --processors 4
    > bioflow analyze --matrix annotome --depth 24 --processors 4

The resulst of analysis will be available in the output folder, and printed out to the standard
output.

For more details or usage as a library, refer to the `usage guide
<http://bioflow.readthedocs.org/en/latest/guide.html#basic-usage>`__.

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
.. |Duplicate Lines| image:: https://img.shields.io/badge/duplicate%20lines-17.66%25-yellowgreen.svg
   :target: http://clonedigger.sourceforge.net/