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
Given the difference in the topology and potential differences in the underlying assumptions, we
pull the interactome knowledge network (where all nodes map to molecular entities and edges - to
physical/chemical interaction between them) and teh annotome knowledge network (where some nodes
might be concepts used to understand the biological systems - such as ontology terms or pathways)
separately.

The parse for interactome is performed by retrieving all the nodes and edges whose ``parse_type``
is ``physical_entity`` for nodes and ``physical_entity_molecular_interaction``, ``identity`` or
``refines``. The giant component of the interactome is then extracted and two graph matrices -
adjacency and laplacian - are build for it. Weights between the nodes are set in an additive
manner according to the policy supplied as the argument to the ``InteractomeInterafce
.full_rebuild`` function or, in a case a more granular approach is needed to the
``InteractomeInterafce.create_val_matrix`` function. By default the
``active_default_<adj/lapl>_weighting_policy`` functions are used from the
``bioflow.algorithms_bank.weigting_policies`` module. Resulting matrices are stored in the
``InteractomeInterface.adjacency_matrix`` and ``InteractomeInterface.laplacian_matrix`` instance
variables, whears the maps between the matrix indexes and maps are stored in the
``.neo4j_id_2_matrix_index`` and ``.matrix_index_2_neo4j_id`` variables.

The parse for the annotome is performed in the same way, but matching ``parse_type`` for nodes to
``physical_entity`` and ``annotation``. In case of a proper graph build, this will result only in
the edges of types ``annotates`` and ``annotation_relationship`` to be pulled. Weighting
functions are used in the similar manner, as well as the mappings storage.


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
In case a specific algorithms to generate pairs of nodes between which
to calculate the information flow is needed, it can be assigned to the ``InteractomeInterface
._flow_calculation_method``. It's call signature should conform to the ``list, list, int ->
list`` signature, where the return list is the list of pairs of ``(node_idx, weight)`` tuples. By
default, the ``general_flow`` method from ``bioflow.algorithms_bank.flow_calculation_methods``
will be used. It will try to match the expected flow calcualtion method based on the parameters
provided (connex within a set if the secondary set is empty/None, star-like if the secondary set
only has one element, biparty if the secondary set has more than one element).

Similarly, methods to evaluate the number of operations and to reduce their number
to a maximum ceiling with the optional int argument ``sparse_rounds`` needs to be assigned to the
``InteractomeInterface._ops_evaluation_method`` and ``InteractomeInterface
._ops_reduction_method``. By default, the are ``evaluate_ops`` and ``reduce_ops`` from
``bioflow.algorithms_bank.flow_calculation_methods``.


Custom random set sampling strategy:
------------------------------------
In case a custom algorithm for the generation of the background sample needs to be implemented,
it should be supplied to the ``InteractomeInterace.randomly_sample`` method as the
``sampling_policy`` argument.

It is expected to accept the an example of sample and secondary sample to match, background from
which to sample, number of samples desired and finally a single string parameter modifying the
way it works (supplied by the ``sampling_policy_options`` parameter of the
``InteractomeInterace.randomly_sample`` method).

By default, this functions implemented by the ``matched_sampling`` fundion in the
``bioflow.algorithms_bank.sampling_policies`` module.


Custom significance evaluation:
-------------------------------
by default, ``auto_analyze`` functions for the interactome and the annotome will use the default
``compare_to_blank`` functions and seek to determine the significance of flow based on comparison
of the flow achieved by nodes of a given degree in the real sample compared to the random "mock"
samples. The comparison will be performed using Gumbel_r function fitted to the highest flows
achieved by the "mock" runs.

As of now, to change the mode of statistical significance evaluation, a user will need to
re-implement the ``compare_to_blank`` functions and mokey-patch them in the modules containing
the ``auto_analyze`` function.
