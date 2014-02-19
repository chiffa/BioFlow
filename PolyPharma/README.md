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

The license is BSD, but in case of academic usage, please cite the url (publication is to come).

Project dependencies:
=====================

**System-level packages:**

* Python 2.7.3 (x86_64 build)
* Java JDK 1.6.0_24 (x86_64 build)
* MongoDB
* Neo4j 1.xx (Note: if you need to install a 2.xx series, please install [gremlin parser engine](https://github.com/neo4j-contrib/gremlin-plugin) too)
* LAPACK, ATLAS, BLAS (only if you are installing Scipy manually)

**Packages installed via Pip:**

* SQLAlchemy 0.8.0b2
* JPype-0.5.4.2 (x86_64 build)
* neo4j-embedded 1.6 (x86_64 build)
* numpy (latest x86_64 build)
* Scipy (latest x86_64 build of LAPACK, ATLAS and BLAS + pip-install)
* bulbs (via pip)
* python-Levenshtein (via pip)
* pymongo
* requests
* scikits.sparse (for the cholesky decomposition of a symetric matrix)
* Sphinx-1.1.3 (documentation build)


Please note that scikits.sparse requires:
 - Cython installation via pip
 - suitesparse-dev(el) installation

This last step is best done via package manager:
On Debian:
```  $ sudo apt-get install suitesparse suitesparse-dev ```
On Fedora / RHCP / CentOS
```  $ sudo yum install suitesparse suitesparse_devel ```

