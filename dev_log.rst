TODOs for the project in the future:
====================================

Confirmed minor refactoring requiring a sane rollback:
------------------------------------------------------

-  Properly indent multi-line :param <parameter type> <parameter name>: descriptors

-  Integrate compops/second estimation to the sources.ini

-  Perform profiling by creating a dedicated set of loggers that would log an "execution time"
flag set

New features:
-------------

-  Add protein abundance level for instantiation of the network

-  Language of network alignment/explanation of net1 by net2: allows
   to compare GO annotation with interactome, cell type specificity 
   analysis or organ context.

-  Add a coarseness feature on the interactome analysis affecting
   sampling behavior, so that precision is sacrificed in favor of
   computation speed.

-  Build an "inspector routine" that would allow us to see the nodes that would allow us to route
   the most information => we need to recompute the most central nodes in Interactome, because we
   still observe a heavy skew in the nodes with a high degree.

-  We always need to first build Interactome Interface before BioKnowledge Interface and in the
   end we need to have both of them build before we can run automated analysis. A nice fix would be
   to raise flags when they are loaded, instead of relying on the loading behaviors.

-  Event sourcing pattern for the graph assembly and modification from the base databases.

-  The execution entry points have to be the five canonical queries.

-  Write the circulation files so that it is possible to calculate the information circulation
   between two sets or as set and a single protein (application for p53 and PKD1 regulators)

-  Distinction between downstream and upstream targets can be implemented by translating the
   directed graph into an associated Markov transition matrix. This will allow to:

    - explicitly allow weight of importance of sources/sinks of information (match the
      distribution shape with the quantile distribution normalization)

    - account only for the information propagating downstream the pathways, not both ways as
      it is the case now. A Markov Matrix differential system solution is a good idea as well, of
      the type F(t) = A*(F(t-1)+B); dA = A*F(t-1)+A*B-F(t)

    - synchronous computation of the flow for all sources/sinks

-  Add citations into the online databases files, that allows integration of different source
   into a single database.

-  Clustering algorithms going beyond the spectral clustering,

    - Not needing a pre-defined number of clusters

    - Able to assign the same node to several clusters

    - Maybe iterative DBSCAN or agglomerative clustering with removal of detected clusters until
      we hit some threshold on the number of number of nodes - average circulation in cluster curve
      obtained from random nodes sampling

    - We can deduce that graph from the clustering of sets of random nodes v.s.

-  Graph exploration module:

    - Strongest eigenvectors / highest circulation in a random set of nodes

    - Randomized clustering/

-  Introduce signal over noise ratio: amount of current in the current
   configuration compared to what we would have expected in case of a
   random set of nodes. We could introduce this as a bootstrap on a
   random subset of nodes to figure which ones are random and which ones
   aren't

-  add @jit() wrapper in order to compile the elements within the current calculation routines.

-  Single command to change the neo4j instance being used or copied

   - Copy a designated database

   - Cd into the designated database and execute neo4j start/stop


Structure-required refactoring:
-------------------------------
-  separate the envelopes for the GO and Reactome graphs retrieval from
   the envelope used to recover and compute over the graph.
   
   -  remove the memoization of individual pairs during the flow withing
      the group computation


-  transfer the annotation search to an ElasticSearch engine.

   -  remove the overhead of loading all the annotation nodes to the
      neo4j instance

   -  allow efficient filtering on the node types. Currently type
      detection and filtering is done upon enumeration. In practice,
      this is not critical, because DB Ids from different databases have
      low intersection

   -  approximate matching capacities for gene names mistypings

-  Split the computation of the blank from the performing of the
   analysis round

-  Inline the background for the InteractomeInstance into the __init__

-  Inline the undumping and dependent variables calculation into the __init__ of
   InteractomeInstance and BioKnowledgeInterface

-  change to element import directories from which too many
   functions/objects are imported (import numpy as np)

-  Make methods running large systems of procedures to being
   dictionary-driven

-  Factor out the traversals used in order to build the Laplacians

-  Refactor the flow calculation as a calculation between two sets protein sets:

    -  Dense calculation or sampling is a strategy

    -  Self-set is just when the two sets are equal

    -  Circulation with a single protein is a special case when one of the sets contain a single
    element.

- Clustering algorithms going beyond the spectral clustering,
    - Not needing a pre-defined number of clusters
    - Able to assign the same node to several clusters
    - Maybe iterative DBSCAN or agglomerative clustering with removal of detected clusters until
    we hit some threshold on the number of number of nodes - average circulation in cluster
    - We can deduce that graph from the clustering of sets of random nodes

Good-to-have; non-critical:
---------------------------

-  In all the DB calls, add a mock-able wrapper that would read the
   state of a project-wide variable and if it is set to True (in
   unittests) will switch it to

-  Bulk-insertion into the neo4j. => Requires taking over the bulbs engine

-  Add active state memoization for the import commander, so that when
   an exception happens, it prints it, terminates gracefully and upon
   restart offers an option to resume from the point of failure while
   managing all the support

-  modify the config generator code so that there is only one place
   where the default configurations are stored and can be modified from
   hte command line interface, instead of a complex CofigsIO class
   management. We actually have several levels of configs:

   -  Configs that are required to properly stitch the code that were
      introduced during the development

   -  Configs managing the third-party services

   -  Configs that are specific to a deploy:

      -  Where to direct the flow of the loggers at every level
      -  Where the datastores are located
      -  How to connect to a database

   -  Configs that allow switching between organisms

      -  Re-filling the database
      -  Re-building the intermediate representations
      -  Re-building the mongoDB reference and average heatmap

-  build a condas-compatible package that would be installable
   cross-plateform and would contain pre-compiled binaries for
   C-extenstions.

-  In case we are calling time-consuming parsers from several location,
   we might want to insert "singleton" module into the block, that
   performs all the parsings only once per program run.

-  We need to dynamically update the values of main_config whenever the location whenever
   configs from configsfiles are modified, so their modification do not require restaring the
   program. Alternatively we say that the configs need to be modified before the rest of the
   program can be executed.

-  Transform all the matrices so that the first one is packed line-based and the second one
   column-based. This will allow the optimization for register pre-pulling in the processor

Possible Major refactoring that would simplify the problem a lot:
-----------------------------------------------------------------

-  Use a dictionary-configurable parser to parse from a given file format to the neo4j database.

    - The dictionary must show what identifiers have to be recovered from the file and to what
        nodes they should be matched in the neo4j database

    - The dictionary must show how the relationships should be inserted into the neo4j database

-  All the insertions are added without node or edge duplication. In case of multiple insertions,
    additional key:value pairs are added to the annotation of the node or the edge

-  Laplacian construction interface takes a dictionary providing instructions on how to
    compute the laplacian or adjacency matrix from the key:map values, both for nodes and edges

    - It allows both easy instantiation from the values for the nodes, such as protein/metabolite
     abundance in a tissue/organ, suppression of an interaction because of a mutation

    - It allows to use a single routine in order to perform different types of computation, such
    as the reliability of information transmission, likehood of randomness/jitter, etc...

    - It allows a high degree of customization by the end user, beyond what would be suggested by
     the initial user

Documentation and description:
==============================

Description:
------------

Following the interaction with Wahid when I was explaining him what my
methods were doing:

-  Explanation of what is current and how uit relations to biology

-  Where are the pathways?

-  Print out the twist ration into the GDF: observed to expected ratio/
   P\_value

- show how to install on a docker and provide the script to perform installation in Ubuntu

- write a quickstart guide

- add pictures of what netowrk analysis looks like

- Validation of results with retrieval of Pamela Silver's paper and John's Overington 300
essential targets: high average information flow and low abundance.

- Generate figures showing the highy-connected nodes in the laplacian matrix corresponding to the
 common chemical molecules (ADP, ATP, Pi, ...). Explain that mechanisms related to such molecules
  would better be described in terms of propositions on actual biological knowledge and that we
  would need to run the two analyses in parallel: both on the concepts and the molecular entities.

- Generate the figures showing that taking in account background that is efficiently reachable by
 a given experimental technique is critical for the proper annotation retrieval, especially for
 the low-informativity terms. Give an example of techniques relying on the aboundance change for
 detection, how they would behave if we randomly sample from back-ground without first setting
 the background.

Internals high-level doc:
-------------------------

-  Limitations: no physical-path toxicity (such as rising pH, changing
   the O2 content or depleting ATP/ADP). They are managed by appropriate GO annotations

-  Retrieving giant connex set and operating on it only.

-  Filtering GOs without enough UP attachement (less than 2) to avoid infinite informativity
   (entropy reduction to 0).


GO Terms analysis techniques
````````````````````````````

-  Perform the statistics on the flow amount and the relation betweeen
   the flow, informativity and confusion potential
-  Perform the statistics on the flow amount and tension for the
   partitions of initial set of proteins to analyse
-  Recover the analysis of the idependent linear groups of the GO terms.
-  Mutual information about the flow and different characteristics, such
   as informativity and confusion potential (which are in fact
   bijective)

Size and memoization pattern of the GO current system:
``````````````````````````````````````````````````````

The current decision is that for the samples of the size of ~ 100
Uniprots, we are better off unpickling from 4 and more by factor 2 and
by factor 10 from 9. Previous experimets have shown that memoization
with pickling incurred no noticeable delay on samples of up to 50 UPs,
but that the storage limit on mongo DB was rapidly exceeded, leading us
to create an allocated dump file.


Specific module improvements:
=============================

This section contains rather general improvements we would like to see in different modules to
make them more independent.

Better data package management:
-------------------------------

Organize the data repository retrieval according to the Python pip convention:

    - use ``package_data`` and ``include_package_data`` to load the pointers to the git
    repositories containing data location.ini files.

    - issue a command to add a git repository mapping a data shortname and data location to a
    downloadable format

    - let user input where the data should be stored on his machine before any actual download
    happens

    - store configuration folders in a ``$HOME/.data_manager/ domain``


Better Reactome parser:
-----------------------

Overall, we want to have a more general and more sane .owl parser

    - Add the parsing of Unification X-REF tags in the Reactome.

    - Unify the parsing structure to the iterative parsing of the tags.

    - Define functions of transformation that will assemble the elements of the owl parsing into
    the class elements. (Flattening the structure)

    - In order to do this, define reduction functions:

        - Inline child's load

        - Discard that attribute

    - The computation of an individual parameter is actually an inlining of a

Beyond something that I am actually needing, this is an excellent exercise at writing a
functional rdf parser that would use a Maybe monad (in case a child/parameter/etc..) is not found

Some of the ideas specific to the bioflow project:

-  perform parsing of unification x-refs in all the meta-types and
   reactions in order to retrieve joins with other databases.

-  return the connecting databases with the number of connections and
   the number of entities getting connected

-  collapse meta-types into a single type and use a type field to
   distinguish them


Better DB_IO management for annot nodes:
----------------------------------------

We want to transfer the load of the indexing to an elasticsearch engine. In order to do that, we
 will suppress the annotation nodes, with their payload and payload typing and transfer it to
 elasticsearch, both with respect to insertion and retrieval. This will allow us to get smaller
 neo4j networks and faster load times.

Beyond that, we would be able to use the mechanism for batch queries on elasticsearch when we are
 retrieving lists to get bulbs identifiers immediately.

Utils and general Utils:
------------------------

**Wrappers:**

-  debug wrapper that logs to the debug channel. In case we are
   performing a graphical debug, we log it as a picture saved with the
   name of calling and the time of calling to the project root

-  visual debugger for the matrix operations that allows to specify what
   input matrixes we would like to inspect and what output matrixes we
   would like to inspect (by index)

Information flow computation:
-----------------------------

**Flow with ponderation**

-  transform the computation to allow for different amount of
   information to be assigned to different nodes.

-  as a rule of thumb, the main computation core does not change, but
   the rules of normalization change.

-  FLOW\_1-2 IMPORTANCE =
   NODE\_1\_IMPORTANCE/TOTAL\_IMPORTANCE\ *NODE\_2\_IMPORTANCE/TOTAL\_IMPORTANCE
   = NODE\_1\_IMPORTANCE*\ NODE\_2\_IMPORTANCE/TOTAL\_IMPORTANCE\*\*2
-  FLOW\_STACK = SUM OF FLOW\_I\_J\_IMPORTANCE\*FLOW(I, J)

**Flow with signs**

-  calculate potentials separately, then perform a summation of
   potentials. Once potentials have been summed, calculate the
   information flow. This however does not reflect much presentation

-  An alternative is to implement a pressure propagation with sign
   inversion to account for positive/negative relations. Even though
   technically relying on the same Laplacian, we will need to
   re-implement routines computing the regulations:

   -  We need to separate reliability flow from the sign propagation
      flow
   -  We would need to enforce the rules that would enforce sign
      propagation only one way: down

-  All in all, we are switching to temperature diffusion on a laplacian
   network. With respect to that, we need a "diffusion" module and a
   separate description of the method how to use it.

**Overall Mathematics**

-  Get rid of Cholesky decomposition: it is not appliable in our
   case because of presence of null eigenvalues In fact there are as
   many null eigenvalues as there are connex segments in the graph

-  Removed: replace pickling by JSON wherever appliable => numpy objects
   are not JSON-seriasable

-  DONE: add the clustering of proteins according to the GO annotation
   similarity

-  TODO: add the evaluations of Zipf-ittude for the proteins

-  DONE: add random matrix filtering-out for the "too noizy" conductions

-  DONE: for the computation of the relevant computational values,
   normalize the connections Graph. Use a laplacian instead of the
   default graph for the decorrelation

-  TODO: add derivatives to analyse scaling factors on for element
   participation in a complex: Is this complex a limiting factor for
   this complex or not?. In case of level variation derivative will be
   the measure for the amount of the trafficked information, whereas in
   case of substantial modification (mutation silencing catalytical
   factor, this will) be the only available one.

-  TODO: add negative/positive potentials for the linkages to the GO
   terms for true Up/Down regulation

-  TODO: orient Zipf-central concepts for different environements (yeah,
   but this is direct biasis, isn't it?) => Better deduce your own
   Zipf-distribution

-  TODO: analyse the sign-connexity of the GO terms analysis tools

-  TODO: add an adaptor for markov model-like analysis - Problem 1: if
   we operate big graphs, we are liklely to run out of memory - Problem
   2: we cannot necessary normalise all the vectors, since some proteins
   are affecting several proteins at the same time

-  TODO: Add the 95% confidence interval for a given percision rate for the depth of sampling.
   For instance if we want 95% confidence into the p_value with 95% confidence, we need to run not
   25 samples, but rather 30 or something in that range.



Features that would be nice to have:
====================================

New analysis features:
----------------------

-  Derivative of GO term flows with respect to a network disruption or protein disruption

-  Negative/positive pressure injection & diffusion in order to account for positive/negative
    regulation in regulatory networks

-  Replace diffusion and flow matrices by causality matrices (directed transitions, allowing to)
    account for upstream/downstream propagation

-  We need to replace eigenmatrix clustering by agglomerative clustring, so that some nodes can
    belong and be important for several clusters instead of having to choose one to which they belong
    more.

-  Stochasticity of transmission: Once we get the abundances of different proteins in the network,


Add protein domain state switches:
----------------------------------

 This will allow us to represent the the changes in protein function following a
 post-translational modification or association in a complex that would be hard to represent
 otherwise.

 More generally, it is switching the distribution of instances between classes that can be
 converted one to another.

Add additional databases:
-------------------------

-  Perform a recovery of post-translational modification sites in
   the normal proteins

-  Perform a recovery of a larger database of the RNAs, both as
   protein transcription elements and as regulatory elements

-  Import the DNA / epigenetic annotation ontology into the
   database to account for the DNA (un)-availability and for the DNA
   transport towards specific (activation or repression regions)

-  Cast in the database Protein Aboundances so that it becomes
   one-and-for-all import Problem: what are we to do in case we are
   willing to use a specific organ and not a general database?

-  Add organ specificity levels of protein expression


Improve crosslinking between different databases
------------------------------------------------

-  Perform a search in the UNIPROT Database in order to improve
   the annotation based on the DisplayNames => this is done separately
   by a matching/lookup module (this would be another good application
   for the elasticsearch engine)

-  Import modification feature from the reactome.org to account for post-translational
    modifications

-  Add fulltext indexes to the nodes (would be another great application of the
   elasticsearch engine)



