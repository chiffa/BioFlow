Data and databases setup:
=========================

Assembling the files required for the database creation:
--------------------------------------------------------

In order to build the main knowledge repository, BioFlow will go and look for the following data
repositories specified in the ``$BIOFLOWHOME/configs/main_configs.yaml`` file:

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

It is possible to specify the file locations and identifiers manually, and then download and
install them. This is to be used when the download locations for the files move.


Similarly, the configs file also controls the organism selection. Three organisms have provided
configurations (human, mouse, S. Cerevisiae). Using the same pattern, other organisms can be
configured, also the lack of data can be a problem (this is already the case for mouse - we
recommend mapping the genes to human if the mouse is used as a model for the organism).


.. WARNING::
    While BioFlow provides an interface to download the databases programmatically, the databases are subject to Licenses and Terms that it's up to the end users to respect


Weighting and secondary set:
============================

In addition to the main analysis set, it is possible to supply to BioFlow as well the secondary
analysis set. In case it is supplied, the information flow will be calculated between the set of
hits and the secondary set. In other terms, only the hypothesis that link the primary set to the
secondary set will be evaluated. This is useful in cases where a set of hits from experimental
results need to be linked to the set of proteins canonically involved in a phenotype or process
of interest.

The secondary set can be supplied to both the ``knowledge_analysis`` and ``interactome_analysis``
as an optional parameter ``secondary_source_list``. This is the ``sec_ids`` return from the
``map_and_save_gene_ids`` method. In order to generate mappints of the ids from the files that
would be compatible with the ``auto_analyze`` methods, there are several possible ways of
supplying file paths::

    > # No background
    > hit_ids, sec_ids, background_ids = map_and_save_gene_ids(hits_file)
    > # Background
    > hit_ids, sec_ids, background_ids = map_and_save_gene_ids(hits_file, background_file)
    > # Single secondary file with background
    > hit_ids, sec_ids, background_ids = map_and_save_gene_ids((hits_file, sec_file), background_file)
    > # mixed list of files with background
    > hit_ids, sec_ids, background_ids = map_and_save_gene_ids([(hits_file, sec_file)
                                                                hits_file_2, ...],
                                                               background_file)

All of them generate valid input for both the ``knowledge_analysis`` and
``interactome_analysis``, but if you are using a list, a list of experiment names is strongly
recommended.


.. WARNING::
    Secondary set has to be disjoing from the primary set. If this is not the case, an explicit error will be raised


For the main as well as the secondary set BioFlow supports weighting. It is possible to supply a



Bioflow supports as well the weighting of the main


Reweighting and node edge exclusion:
====================================