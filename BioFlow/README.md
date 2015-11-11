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

The license is BSD, but in case of academic usage, please cite the *url* (publication is to come).



Project dependencies:
=====================

**System-level packages:**

* Python 2.7.x (x86_64 )
* Java JDK 1.6/1.7
* Neo4j 1.xx (Note: if you need to install a 2.xx series, please install [gremlin parser engine](https://github.com/neo4j-contrib/gremlin-plugin) too)
* MongoDB
* suitesparese & suitesparse-dev
* LAPACK, ATLAS, BLAS (only if you are installing Scipy manually)


**Packages installed via Pip:**

* JPype (x86_64 build)
* Cython
* neo4j-embedded (x86_64 build)
* SQLAlchemy (latests x86_64 build)
* NumPy (latest x86_64)
* SciPy (latest x86_64)
* bulbs (via pip)
* python-Levenshtein (via pip)
* pymongo
* requests
* Sphinx (documentation build)
* scikits.sparse (for the cholesky decomposition of a symetric matrix)


Please note that scikits.sparse requires Cython installation via pip and suitesparse-dev(el) installation.

This last step is best done via package manager:

On *Debian*:   ```  $ sudo apt-get install suitesparse suitesparse-dev ```

On *Fedora* / *RHCP* / *CentOS*:    ```  $ sudo yum install suitesparse suitesparse_devel ```


Required files:
===============
I will write a script to download and install those files automatically from within the application later on,
provided that there are some renaming conventions


**Databases:**
* OBO 1.2 file of GO terms and relations, available at: http://www.geneontology.org/GO.downloads.ontology.shtml
* UNIUPROT-SWISSPROT text database file available at: http://www.uniprot.org/downloads
* Reactome.org "Events in the BioPax level 3" file, available at: http://www.reactome.org/download/index.html
* HiNT binary interaction files for the organisms of interest, availble at: http://hint.yulab.org/batch.html
* BioGRID files in the tab2 format, available at http://thebiogrid.org/download.php
* Gene-chromosome mapping files from the Uniprot documentation: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/
* Organism-specific protein aboundance files, available at: http://pax-db.org/#!downloads

Please, pay attention to correctly rename the files or edit the files names so that the they correspond to the sources.ini
file in the /configs directory of the application
