PolyPharma
=======================

This code is currently under development and thus it isn't always stable, thoroughly tested or well documented.

This project's goal is to predict systemic effect of multiple gene perturbation, wherther trigered by a drug or by
a disease (such as cancer or disease with complex genetic background). It's main intended uses are prediction of
drug toxicity of de-novo drugs due to a distributed off-target effect and linkage between a phenotype and a complex
genotype.

It's main advantage is integration of quantitative computational predictions with prior biological knowledge and
ability to integrate such diverse source of knowledge as databases, simulation, publication data (currently in dev)
and expert knowledge.

This application is provided as-is, without any warranties or support. Use it at your onw risk

However, if you are willing to test it and encounter problems or are willing to provide feedback, please fill in
an issue ticket on GitHub and I will be glad to assist you in the measure of my possibilities.

Project dependencies:

Python 2.7.3 (x86_64 build)
SQLAlchemy 0.8.0b2
JPype-0.5.4.2 (x86_64 build)
Java JDK 1.6.0_24 (x86_64 build)
neo4j-embedded 1.6 (x86_64 build)
numpy (latest x86_64 build)
Scipy (latest x86_64 build of LAPACK, ATLAS and BLAS + pip-install)
bulbs (via pip)
python-Levenshtein (via pip)
scikits.sparse (for the cholesky decomposition of a symetric matrix)
pymongo (and mongodb database running on the default port)

Please note that scikits.sparse requires:
 - Cython installation via pip
 - suitesparse-dev install and CHOLMOD .so files linked and referenced by the compiler.

 Currently the only reasonably easy way to get scikits to compile agains suitesparse is on Debian Linux and
 relies on issuing an:
  $ sudo apt-get install suitesparse suitesparse-dev


Documentation Build:
Sphinx-1.1.3