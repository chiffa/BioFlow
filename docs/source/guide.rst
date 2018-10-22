Data and databases setup:
=========================

Assembling the files required for the database creation:
--------------------------------------------------------

In order to build the database, the program is going to look for the following files specified
in the following locations within the bioflow/configs/sources.ini:

* OBO 1.2 file of GO terms and relations,

    * dowloaded from: `here <http://purl.obolibrary.org/obo/go/go-basic.obo>`__
    * will download/look for for go.obo file at $DB_HOME$/GO/

* UNIPROT-SWISSPROT .txt text database file

    * downloaded from `here <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz>`__
    * will store/look for the master tab file at $DB_HOME$/Uniprot/uniprot_sprot.dat
    * will load the information specified by the NCBI tax id for the organism currently loaded

* Reactome.org "Events in the BioPax level 3" file

    * downloaded from `here <http://www.reactome.org/download/index.html>`__
    * will store/look for the files relevant to the organism at $DB_HOME$/Reactome/

* HiNT PPI binary interaction files for the organisms of interest

    * dowloaded from `here <http://hint.yulab.org/download/HomoSapiens/binary/hq/>`__, `here <http://hint.yulab.org/download/SaccharomycesCerevisiaeS288C/binary/hq/>`__ and `here <http://hint.yulab.org/download/MusMusculus/binary/hq/>`__
    * will store/look for the files SaccharomycesCerevisiaeS288C_binary_hq.txt, HomoSapiens_binary_hq.txt, MusMusculus_binary_hq.txt at the path $DB_HOME$/HiNT/

* BioGRID ALL_ORGANISMS PPI file in the tab2 format

    * dowloaded from `here <http://thebiogrid.org/download.php'>`__
    * will store/look for the files at $DB_HOME$/BioGRID/

* TRRUST literature-curated TF-target interaction files in tsv format

    * downloaded from `here <http://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv>`__ and `here <http://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv>`__
    * will store/look for the files at $DB_HOME$/TFs/TRRUST

* IntAct ComplexPortal tsv files

    * dowloaded from `here <ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/homo_sapiens.tsv>`__ and `here <ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/saccharomyces_cerevisiae.tsv>`__
    * will store/look for the files at $DB_HOME$/ComplexPortal

* Phosphosite protein kinase-subrstrate tsv files

    * downloaded from `here <https://www.phosphosite.org/staticDownloads>`__
    * will store/look for the files at $DB_HOME$/PhosphoSite

It is possible to specify the file locations and identifiers manually, and then download and install them
to the specified locations manually.

However the following command should be able to do it for you for three commonly used organism (human, mouse, saccharomyces cerevisae),
provided you follow the instructions properly::

    > python CLUI.py initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver mongodb://localhost:27017/

    > python CLUI.py downloaddbs

    > python CLUI.py setorgconfs --organism [mouse, human, yeast]

.. WARNING::
    While BioFlow provides an interface to download the databases programmatically, the databases are subject to Licenses and Terms that it's up to the end users to respect


Typical sources.ini configfile:
-------------------------------

Here is what a typical configfile would look like::

    [HINT]
    file = HomoSapiens_binary_hq.txt

    [BIOGRID]
    name_pattern = Homo_sapiens

    [INTERNAL]
    mongosuffix = _v_1
    dumpprefix = /human
    mongoprefix = _human
    compops = 1

    [UNIPROT]
    tax_ids = 9606

    [GO]
    file = go.obo

    [REACTOME]
    file = Homo_sapiens.owl

    [TRRUST]
    file = trrust_rawdata.human.tsv
    significance = 1

    [COMPLEXPORTAL]
    file = homo_sapiens.tsv

    [PHOSPHOSITE]
    file = Kinase_Substrate_Dataset.tsv
    organism = human



The configuration files might be declared and switched manually (only the "source.ini" one will be parsed,
folders such as "sources_organism.ini" will be ignored and can be renamed to "source.ini" quite easily)

It is possible for the users to generate source.ini file for three organisms with the following command::

    python CLUI.py setorgconfs --organism [mouse, human, yeast]

This allows to switch rapidly between different investigated organism.

Please don't forget to switch or purge neo4j databases between organisms, because each organism needs it's own neo4j instance.