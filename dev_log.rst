This file regroups the TODOs for the project in the future.
===========================================================

Confirmed minor refactoring requiring a sane rollback:
------------------------------------------------------

-  render the usage of analytical/background sets uniform throughout the
   system.

-  split the computation of the blank from the performing of the
   analysis round

Minor refactoring suggestions:
------------------------------

-  add active state memoization for the import commander, so that when
   an exception happens, it prints it, terminates gracefully and upon
   restart offers an option to resume from the point of failure while
   managing all the support

-  separate configs by the level of modification frequency and move to
   the top-level directory

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

      -  Forking and switching the databases
      -  Re-filling the database
      -  Re-building the intermediate representations
      -  Re-building the mongoDB reference and average heatmap

-  build a condas-compatible package that would be installable
   cross-plateform and would contain pre-compiled binaries for
   C-extenstions.

-  In case we are calling time-consuming parsers from several location,
   we might want to insert "singleton" module into the block, that
   performs all the parsings only once per program run.

-  In case of go/uniprot and reactome parsers, create wrapping objects,
   that would contain memoization dictionaries.

-  change to element import directories from which too many
   functions/objects are imported (import numpy as np)

-  refactor the configs to contain the following elements:

   -  Configs computed based on the .ini files
   -  Inner configs that has an effect on what is being computed
   -  Inner configs that has only anything to do with the project
      structure and implementation

-  remove the memoization of individual pairs during the flow withing
   the group computation

-  Make methods running large systems of procedures to being
   dictionary-driven

-  create a set of wrappers that would log, either into INFO or into
   DEBUG channel, every time function is called and parameters on the
   entrance and the exit

-  

Major refactoring suggestions:
------------------------------

-  transfer the annotation search to an ElasticSearch engine. Reasons:

   -  remove the overhead of loading all the annotation nodes to the
      neo4j instance

   -  allow efficient filtering on the node types. Currently type
      detection and filtering is done upon enumeration. In practice,
      this is not critical, because DB Ids from different databases have
      low intersection

   -  approximate matching capacities for gene names mistypings

-  in all the DB calls, add a mock-able wrapper that would read the
   state of a project-wide variable and if it is set to True (in
   unittests) will switch it to

-  Introduce signal over noise ratio: amount of current in the current
   configuration compared to what we would have expected in case of a
   random set of nodes. We could introduce this as a bootstrap on a
   random subset of nodes to figure which ones are random and which ones
   aren't

-  Add protein abundance level

-  Add a coarseness feature on the interactome analysis affecting
   sampling behavior, so that precision is sacrificed in favor of
   computation speed.

-  separate the envelopes for the GO and Reactome graphs retrieval from
   the envelope used to recover and compute over the graph.

-  

Required documentation:
-----------------------

Following the interaction with Wahid when I was explaining him what my
methods were doing:

-  Explanation of what is current and how uit relations to biology

-  Where are the pathways?

-  Print out the twist ration into the GDF: observed to expected ratio/
   P\_value

- show how to install on a docker and provide the script to perform installation in Ubuntu

- write a quickstart guide

- add pictures of what netowrk analysis looks in

Improvements to different modules:
----------------------------------

Biological data parsers
~~~~~~~~~~~~~~~~~~~~~~~

**Reactome BioPax lvl3 owl parser:**

-  perform parsing of unification x-refs in all the meta-types and
   reactions in order to retrieve joins with other databases.

-  return the connecting databases with the number of connections and
   the number of entities getting connected

-  collapse meta-types into a single type and use a type field to
   distinguish them

Utils and general Utils:
~~~~~~~~~~~~~~~~~~~~~~~~

**Wrappers:**

-  debug wrapper that logs to the debug channel. In case we are
   performing a graphical debug, we log it as a picture saved with the
   name of calling and the time of calling to the project root

-  visual debugger for the matrix operations that allows to specify what
   input matrixes we would like to inspect and what output matrixes we
   would like to inspect (by index)

Information flow computation:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Add additional Sources/Dimensions
---------------------------------

-  DONE: perform a recovery of important domains from PDB

-  TODO: perform a recovery of post-translational modification sites in
   the normal proteins

-  TODO: Perform a recovery of a larger database of the RNAs, both as
   protein transcription elements and as regulatory elements

-  TODO: Import the DNA / epigenetic annotation ontology into the
   database to account for the DNA (un)-availability and for the DNA
   transport towards specific (activation or repression regions)

-  TODO: Cast in the database Protein Aboundances so that it becomes
   one-and-for-all import Problem: what are we to do in case we are
   willing to use a specific organ and not a general database?

-  DONE: Add ENSEMBL idnetifiers and gene names indexing

-  DONE: In the Uniprot insertions, switch from the hard filtering
   (inserting only uniprots with acnums accessible) from the reactome to
   ALL the uniprots, but using the "inclusion parameter.

-  TODO: add organ specificity levels of protein expression

-  Rejected: remove hard filtration on too participative nodes; instead
   treat it with variational coefficients => Excessively increases
   complexity

Improve crosslinking between different databases
------------------------------------------------

-  TODO: perform a search in the UNIPROT Database in order to imoprove
   the annotation based on the DisplayNames => this is done separately
   by a matching/lookup module

-  TODO: we might want to parse the traceability of the all the
   compounds and link by adding the xref parsed information to them.
   This might be critical to adress the issues imposed by the difference
   in the database versions

-  TODO: import modification feature insertion from the reactome
   database to account for post-translational modifications

-  TODO: verify if GO\_Terms analysis conserves the "regulation
   relations or not.

-  TODO: add fulltext indexes to the nodes

-  TODO: There might be an error in the module responsible for linkage
   between the uniprots and the accession numbers: for instance the
   20253 has an annotation with an Acnum, but has no Uniprot attached to
   it within the database => this is possibly due to the fact that some
   of the uniprots are refered as being from different organsims (such
   as HIV invasion pathway)

From the mathematical point of view
===================================

-  TODO: Get rid of Cholesky decomposition: it is not appliable in our
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

Utils module
------------

-  TODO: In the xml\_doctor, add correlation between presencesof
   different subtypes in the references?

Current Limitations
-------------------

-  Limitations: no physical-path toxicity (such as rising pH, changing
   the O2 content or depleting ATP/ADP)

Potential enhancements:
-----------------------

-  TODO: create GO and Pathway Structure access

   -  Calibrate the values so that after ~ 3 transitions the correlation
      vanishes on average (Follow Pamela Silver Approach) => this is
      actually the cumulated perturbation of
   -  two targets that should vanish totally

-  TODO: along with Overingtonicity integrate the list of essential
   genes in human diseases from the PLoS 2011 publication

   -  Importance of complementation of the information with the
      Reactome.org data with the EHiT data: otherwise the information
      circulation completely sucks
   -  Reactome.org: the interactions due to kinases aren't explicitly
      shown. Instead a broadcasting through the secondary features that
      perform the modification
   -  Is needed. Which is completely stupid, because it doesn't show the
      specific action on the proteins due to the conformation
      modification. Thus Reactome.org
   -  is more of a ressource for human experts then for truly
      machine-learning tasks.

To be treated:
--------------

::

    # If a specific set of GO_Terms is put down, we can say that the function they describe is down.
    # Recall v.s. precision for a GO array for a perturbed protein set?
    # Non-randomness of a recall?
    # Pathway structure?

    # Method extendable to inhibition / activation binaries, by introducing positive / negative values for the matrix

    # Fill in the matrix with the values
    # Take an impact vector
    # Continue multiplications as long as needed for convergence

    # export the matrix as a flat file
    #    => Most significantly touched elements, especially in the UNIPORT
    #    => Get the vector of affected proteins, then multiply it over the transfer
    #        Matrix until an equilibrium is reached.

    # Pay attention to the criticality spread => vector shoud increase exponentially for the important prots, effectively shutting down the whole system
    # But not in the case of "unimportant proteins"

    # => Assymetric influence matrices (causality followship)
    # Markov clustering linalgebra on sparce matrices to accelerate all this shit?

    # We could actually envision it as a chain reaction in a nuclear reactor, leading either to a reaction spiraling out of control (total functional shutdown, at least for a
    # given function.

    # Idea behind the eigenvectors: if we generate random sets of genes perturbating the network, some combination would lead to a way more powerful effect when propagated
    # in a markovian, turn-based network (runaway), whereas other sets will lead to a lighter runaway. A way to estimate runaway specifics of protein-protein interaction network
    # The strongest runaway would be generated by the highest absolute-value link

    # Group node definintion have to be corrected so they are not all related together but instead are linked towards the central "group" node!!!!


    # Shut down HiNT analysis => Slightly improves the result

    # Synchronious eigenvectors approach: protect agains entering into a forbidden list the target node
    # start iterating matrix multiplications starting from the node1 to go to the node2
    # enter each node visited in the forbidden set, except for node2
    # terminate iterating when there are no more new reaches for node2 after all the interations

    # Percentage of information reaching a given node compared to all the information reaching the node: eigenvalue approach too.
    # Error we do: compute three times

    # Ok, what is going on is that we have collections of ~ 300 elements completely screwing our system

    # The problem that a information broadcasting between the elements of the same group is not a good thing, but a direct broadcasting into a reaction is actually
    # what we need in our matrix.


    # In order to be precise, we should not only take in account the power of bindinb between a molecule and protein and criticality of the protein, but also the abundance of the
    # protein in the reactome

    # => Done with the aboundance retrieval

    # DONE: use sparse matrixes routines to calculate the number of connex elements in the graph
    #   Problem: there are 58 disconnected sets.
    #   Solution: retrieve the Node Ids of the main connex Set and write them into the neo4j graph, then retrieve only them

    # DONE: markup of the major connex graph within neo4j database
    #    Waiting for the execution


    # DONE: calculate the distance graph
        # seems to work pretty well with Djikistra.
        # Can we perform a retrieval of specific nodes within distance X of the main component?

    # DONE: buid jump tables to compute the number of reactional transitions
    #    Implemented by using djikstra algo from scipy.sparse.csgraph
    #

    # DONE: retrieve Pamela silver's degradation of the data with the time
    #    Waiting for the execution
    #

    # DONE: pull in the annotations regarding the proteins aboundances
    #
    #

    # DONE: pull in the 300 essential targets from the EBI dude (John Overington)
    #     Results aren't so conclusive. It seems that the protein concentration defenitely plays some role in the determining if a protein is a
    #     Target of an existing drug or not, butthe informativity seems not. Probably this is due to the fact that the targeted proteins are often
    #     cellular receptors.

    # DONE: perform a localization factor pull-out for the Uniprots based on their proteins of attachement
    #        Waiting for the execution

    # DONE: broadcast to uniprots for the localization of the pointed proteins

    # # DONE: reverse GO_Access: provided the Uniprots find the proteins carrying over the most information
    # DONE: mount a PyMongo data store in order to be able to save and retrieve the programming objects easily
    #         How is it done: - picket to string
    #        Store an object in a collection defined by it's Id and computation number
    #        If requested, retrieve by ID or else
    #         Index on the GO ID and belonging UNIPROTs (If same set of uniprots, it is the same) => store as sets
    #         Pickles of sets with the same elements are always the same

    # DONE: remake the sampling so it is efficiently 170**2/2 one to one randomly chosen pairs that are calculated, and not the whole 170 ensemble, so that the
    # Informativities actually follow a gaussian distribution

    # DONE: filter out GOs with not enough UP

    # DONE: export of the analytical system in a Gephy-compatible GDF
    #       => Yes, export as GDF, including attached proteins, with names and GOs with informativities and random pick orobas

General programming:
--------------------

Unit-testing
~~~~~~~~~~~~

-  Create a whack xml, then run all the database loads/unloads one after
   another to check if everything is present and is working as expected.
-  Create smaller unit-tests to check if matrix manipulations work
   correctly

Traceback of programming decisions:
-----------------------------------

GO Analysis and visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GO Terms analysis techniques
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Perform the statistics on the flow amount and the relation betweeen
   the flow, informativity and confusion potential
-  Perform the statistics on the flow amount and tension for the
   partitions of initial set of proteins to analyse
-  Recover the analysis of the idependent linear groups of the GO terms.
-  Mutual information about the flow and different characteristics, such
   as informativity and confusion potential (which are in fact
   bijective)

Size and memoization pattern of the GO current system:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current decision is that for the samples of the size of ~ 100
Uniprots, we are better off unpickling from 4 and more by factor 2 and
by factor 10 from 9. Previous experimets have shown that memoization
with pickling incurred no noticeable delay on samples of up to 50 UPs,
but that the storage limit on mongo DB was rapidly exceeded, leading us
to create an allocated dump file.
