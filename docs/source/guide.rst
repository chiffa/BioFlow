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

Adding new data to the main knowledge repository:
-------------------------------------------------
The easiest way to add new information to the main knowledge repository is by finding the nodes
to which new knowledge will attach (provided by the ``convert_to_internal_ids`` function from the
``bioflow.neo4j_db.db_io_routines`` module for a lot of x-ref identifiers for physical entity
nodes), and then process to add new relationships and nodes
using the functions ``DatabaseGraph.link`` to add a connection between nodes and ``DatabaseGraph
.create`` to add a new node. ``DatabaseGraph.attach_annotation_tag`` can be used in order to
attach annotation tags to new nodes that can be searcheable from the outside. All functions can
be batched (cf API documentation).

A new link will have a call signature of type ``link(node_id, node_id, link_type, {param: val})
``, where node_ids are internal database ids for the nodes provided by the
``convert_to_internal_ids`` function, ``link_type`` is a link type that would be handy for you to
remember (preferably in the snake_case). Two parameters are expected: ``source`` and
``parse_type``.  ``parse_type`` can only take a value in ``['physical_entity_molecular_interaction',
'identity', 'refines', 'annotates', 'annotation_relationship', 'xref']``, with ``'xref'`` being
reserved for the annotation linking.

A new node will have a call signature of type ``create(node_type, {param:val})`` and return the
internal id of the created node. ``node_type`` is a node type that would be handy for you to
remember (preferably in the snake_case). Four paramters are expected: ``'parse_type'``,
``'source'``, ``'legacyID'`` and ``'displayName'``. ``'parse_type'`` can take only values in
``['physical_entity', 'annotation', 'xref']``, with ``'xref'`` being reserved for the annotation
linking. ``legacyID`` is the identifier of the node in the source database and ``displayName`` is
the name of the biological knowledge node that that will be shown to the end user.


Main knowledge graph parsing:
-----------------------------



Custom weighting function:
--------------------------
In order to account for different possible considerations when deciding which nodes and
connections are more likely to be included in hypothesis generation, we provide a possibility for
the end user to use their own weight functions for the interactome and the annotome.

The provided functions are stored in ``bioflow.algorithms_bank.weighting_policies`` module. An
expected signature of the function is ``starting_node, ending_node, edge > float``, where
``starting_node`` and ``ending_node`` are of ``<neo4j-driver>.Node`` type, whereas ``edge`` is of
the ``<neo4j-driver>.Edge`` type. Any properties available stored in the main knowledge
repository (neo4j database) will be available as dict-like properties of node/edge objects
(``<starting/ending>_node['<property>']``/``edge['property']``).

The functions are to be provided to the ``bioflow.molecular_network
.InteractomeInterface.InteractomeInterface.create_val_matrix()`` method as
``<adj/lapl>_weight_policy_function`` for the adjacency and laplacian matrices respectively.


Custom flow calculation function:
---------------------------------


Custom random set sampling strategy:
------------------------------------


Custom significance evaluation:
-------------------------------

