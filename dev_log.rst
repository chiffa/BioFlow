TODOs for the project in the future:
====================================

On the table:
-------------
<Documentation>

 - DONE: [DOC] pass and APIdoc all the functions and modules

 - DONE: [DOC] document all the possible exceptions that can be raised
    - what will raise an exception

 - DONE: [SANITY] remove all old dangling variables and code (deprecated X)

 - TODO: [DOC] Document the proper boot cycle of the application
    - $BIOFLOWHOME check, use the default location (~/bioflow)
    - in case the user configs .yaml is not found, copy it from its own registry to $BIOFLOWHOME
    - use the $BIOFLOWHOME/config/main_config.yaml to populate variables in main_configs
    - use the information there to load the databases from the internet and set the parsing
        locations
    - be careful with edits - configs are read safely but are not checked, so you can get random
        deep python errors that will need to be debugged.
    - everything is logged to the run directory and the $BIOFLOWHOME/.internal/logs

 - TODO: [DOC] put an explanation of overall workflow of the library

 - DONE: [DOC] check that all the functions and modules are properly documented

 - DONE: [SANITY] move the additional from the "annotation_network" to somewhere saner =>
        Separate application importing BioFlow as a library

 - TODO: [DOC] document where the user-mapped folders live from Docker

 - TODO: [DOC] document the user how to install and map to a local neo4j database

 - TODO: [REFACTOR] re-align the command line interface onto the example of an analysis pipeline

 - TODO: [SANITY] Docker:
    - Add outputs folder map to the host filesystem to the docker-compose
    - Remove ank as point of storage for miniconda in Docker



Current refactoring:
--------------------

 - TODO: [KNOWN BUG] [CRITICAL]: pool fail to restart on cholmod usage
    There seems to be an interference between the "multiprocessing" pool method and a third-party
    library (specifically cholmod). It looks like after a first pool spawn, cholmod is failing to
    re-load a new solution again.
        - DONE: try re-loading the problematic module every time => Fails
        - DONE: try to explicitely terminate the pool. It didn't work in the end
        - TODO: extract a minimal example and post it to StackOverflow
        - DONE: add an implicit switch if there is a single element to analysis or multiple
            between explicitely multi-threaded and implicitely single-threaded
        - DONE: temporary patch and flag it as a known bug

 - TODO: [FEATURE]: Factor out the structural analysis of the network properties to a module

 - TODO: [FEATURE]: Factor out the clustering analysis of the network to a different function in
        the knowledge/interactome analyses

 - TODO: [REFACTOR]: factor out the process spawning logic shared between knowledge and
        interactome analysis to a "utility" domain of BioFlow

 - TODO: [USABILITY]: adjust the sampling spin-up according to how many "good" samples are
        already in the mongodb

 - TODO: [SANITY] [REFACTOR]: sanify the databases management
    - TODO: create the `data_stores` package and move everything relating to mongoDB and neo4j there
    - TODO: remove the `GraphDeclarator.py`, fold the logic directly into the `cypher_drivers`
    - TODO: rename `db_io_routines` to `bionet_routines`
    - TODO: import `cypher_drivers` as `bionet_store`
    - TODO: rename `cypher_drivers` to `_neo4j_backend`
    - TODO: import the `mongodb.py` as an alias with `samples_storage`
    - TODO: fold the laplacians .dump object storage in dumps as `auxilary_data_storage`
    - TODO: put a type straight-jacket
    - TODO: move the `internet_io` to the `data_stores` package

 - TODO: [SANITY]: put all the imports  under the umbrella making clear where they come from

<Sanify BioKnowledge>

 - TODO: Develop a pluggable Informativity weighting function for the matrix assembly

 - TODO: Allow for a score for a physical entity term attachment to the ontology system
    - eg. GO attachment comes from UNIPROT
    - reactome comes from reactome and can be assigned a linkage score.

 - TODO: inline the reach computation to remove excessively complex function

<Weighting of the nodes>:

 - TODO: Define pairs in the sampling with a "charge" parameter if the parameters supplied by the
    - The problem is that there is no good rule for performing a weight sampling, given that there
        are now two distribution in the interplay
    - We however cannot ignore the problem, because we discretize a continuous distribution -
        something that is a VERY BAD PRACTICE (TM)
    - Basically, the problem is how to perform statistical tests to make sure not to make
        overconfident calls.
        => degree vs weight - based sampling?

 - TODO: allow for weighted sets and biparty sets to be computed:
    - TODO: modify the flow computation functions to allow for a current to be set
    - TODO: allow the current computation to happen on biparty sets and weighted sets
    - TODO: perform automated switch between current computation policies based on what is
        supplied to the method

 - TODO: allow for different background sampling processes
    - TODO: add two additional parameters to the mongoDB:
        - TODO: parameter of the sampling type (set, weighted set, biparty, weighted biparty)
        - TODO: parameter specific to the set type: (set size, set size + weight distribution,
                                                    pairs, pairs + weights)
    - TODO: add a sampling policy transformer that takes in the arguments and returns a proper
        policy:
        - DONE: set sampling
        - TODO: weighted set sampling
        - TODO: bipary sampling (~ set sampling)
        - TODO: weighted biparty sampling (~ weighted set sampling)


<Type hinting and typing>
 - TODO: [MAINTENABILITY] put a straightjacket of the types of the tuples passed around and
        function type signatures => Partially done, long-term project

 - TODO: [SANITY] convert the dicts to Type Aliases / NewType and perform proper type hinting
        => Partially done, long-term project

 - TODO: [SANITY]: define appropriate types: => Partially done, long-term project
    - neo4j IDs
    - laplacian Matrix
    - current
    - potential


<Environment registration>

 - TODO: [USABILITY] store a header of what was analyzed and where it was pulled from + env
        parameters in a text file in the beginning of a run.

 - TODO: define a persistent "environment_state" file in the $BIOFLOWHOME/.interal/
    - TODO: log the organism currently operating
        => check_organism() -> base, neo4j, laplacians
        => update_organism(base, neo4j, laplacians) -> None
    - TODO: log the organism loaded in the neo4j
    - TODO: log the organism loaded in laplacians
    - TODO: define a "check_org_state" function
    - TODO: define an "update_org_state" function
    - TODO: make sure that the organisms in operating/neo4j/laplacian are all synced
        => "check_sync()" (calls check_organism, raises if inconsistency)
    - TODO: make sure that the neo4j is erased before a new organism is loaded into it.
        => "check_neo4j_empty" (calls check_organism, checks that neo4j is "None")
    - TODO: make sure that the retrieved backtround set is still valid

<Memoization of actual analysis runs>

 - TODO: [FEATURE] currently, re-doing an analysis with an already analyzed set requires a complete
        re-computation of flow generated by the set. However, if we start saving the results into
        the mongoDB, we can just retrieve them, if the environement and the starting set are
        identical.

 - TODO: [USABILITY] move the dumps into a mongo database instance to allow swaps between builds
        - wrt backgrounds and the neo4j states

 - TODO: [USABILITY] since 4.0 neo4j allows multi-database support that can be used in order to
    build organism-specific databases and then switch between them, without a need to rebuild

 - TODO: [USABILITY]: allow a fast analysis re-run by storing actual UP groups analysis in a
        mongo database - aka a true memoization.


 - ????: [USABILITY] add the Laplacian nonzero elements to the shape one (????)


<pretty progress>

 - TODO: [USABILITY] Improve the progress reporting
        Move the INFO to a progress bar. The problem is that we are working with multiple threads in
        async environment. This can be mitigated by using the `aptbar` library
    - TODO: single sample loop to aptbar progress monitoring
    - TODO: outer loop (X samples) to aptbar progress monitoring
    - TODO: move parameters that are currently being printed in the main loop in INFO channel to
        DEBUG channel
    - TODO: provide progress bar binding for the importers as well

 - TODO: [USABILITY]: fold the current verbose state into a `-v/--vebose` argument


 - DONE: perform the explicit background pass for BioKnowledge as for the Interactome

 - DONE: rename the 'meta_objects' to 'Reactome_base_object' in reactome_inserter.py

<node weights/context forwarding>

 - DONE: eliminate InitSet saving in KnowledgeInterface and Interactome interface
    - DONE: they become _background
    - DONE: _background can be set on init and is intersected with accessible nodes (that are saved)
        - on _init
        - on _set_sampling_background

 - DONE: perform the modification of the background selection and registration logic.
        DONE: First, it doesn't have to be integrated between the AnnotationAccess interface and
            ReactomeInterface
        DONE: Second, we can project the background into what can be sampled instead of
            re-defining the root of the sampling altogether.
        DONE: Background is no more a parameter supplied upon constructions, but only for the
            sampling, where it still gets saved with the sampling code.
        DONE: the transformation within the sampling is done

 - DONE: deal with the parse_type inconsistencies (likely remainders of the previous insertions
that was off)
    - DONE: (physical_entity)-refines-(annotation)
    - DONE: (annotation)-refines-(annotation)
    - DONE: (annotation)-annotates-(annotation)
    - DONE: is_next_in_pathway still has custom_from and custom_to

 - DONE: change the way the connections between GOs and UPs are loaded into the KnowledgeInterface

 - DONE: [REFACTOR] The policy for the building of a laplacian relies on neo4j crawl (2 steps)
    and the matrix build:
    - neo4j crawl
        - A rule/routine to retrieve the seeds of the expansion
        - A rule/routine to expand from those routines and insert nodes into the network
    - matrix build
        - creates the maps for the names, ids, legacy IDs and matrix indexes for the physical
            entities that will be in the interactome
        - connect the nodes with the links according to a weighting scheme
        - normalize the weights for the laplacian

 - DONE: STAGE 2/3:
    - DONE: parse the entire physical entity graph
    - DONE: convert the graph into a laplacian and an adjacency matrix
    - DONE: check for the giant component
    - DONE: write the giant component
    - DONE: re-parse the giant component only
    - DONE: re-convert the graph into a laplacian

 - TODO: ax the deprecated code and class variables in InteractomeInterface

 - TODO: ax the deprecated variables in the internal_configs


 - DONE: [REFACTOR] inline the neo4j classes deletion(the same way as self_diag)

 - DONE: [REFACTOR] On writing into the neo4j DB we need to separate the node types and edges:
    - node: physical entity nodes
    - edge: physical entity molecular interaction
    - edge: identity
    - node: annotation (GO + Reactome pathway + Reactome reaction)
    - edge: annotates
    - node: x-ref (currently the 'annotation' node)
    - edge: reference (currently the 'annotates' edge type)
    For compatibility with life code, those will initially be referred to as parse_types as a
        property

 - DONE: [REFACTOR] add universal properties:
    - N+E: parse_type:
        - N:
            - physical_entity
            - annotation
            - xref
        - E:
            - physical_entity_molecular_interaction
            - annotates
            - annotation_relation
            - identity
            - reference
            - refines
    - N+E: source
    - N+E: source_<property> (optional)
    - N: legacyID
    - N: displayName

 - DONE: [REFACTOR] check that the universal properties were added with an exception in
    - DB.link if `parse_type` not defined or `source` not defined
    - DB.create if `parse_type` not defined, `source` not defined, `legacyID` not defined or
        `displayName` not defined

 - DONE: [REFACTOR]
    - NOPE: Either add a routine that performs weight assignment to the nodes
    - DONE: Or crawl the nodes according to the parse_type tags, return a dict of nodes and a
        dict of relationships of the types:
            - NodeID > neo4j.Node
            - NodeID > [(NodeID, OtherNodeID), ] + {(NodeID, OtherNodeID): properties}


 - DONE: [DEBUGGING] write a tool that checks the nodes and edges for properties and numbers and
    then finds patterns
        - DONE: nodes
        - DONE: edges
        - DONE: patterns
        - DONE: formatting

 - DONE: PROBLEM:
    - 'Collections' are implicated in reactions, not necessarily proteins themselves.
    - Patch:
        - Either: link the 'part of collection' to all the 'molecular entity nodes'
        - Or: create 'abstract_interface'
        - Or: same
    - DONE: Due to a number of inclusions (Collection part of Collection, ....), we are going to
        introduce a "parse_type: refines"

Current rewriting logic would involve:
    - DONE: Upon external insertion, insert as well the properties that might influence the
    weight computation for the laplacian construction
        - cross-link the Reactome nodes linked with a "reaction" so that it's a direct linking in
        the database
    - DONE: Changing neo4j crawl so that it uses the edges properties rather than node types
        - For now we will be proceeding with the "class" node properties as a filter
        - Crawl allowed to pass through edges with a set of qualitative properties
        - Crawl allowed to pass through nodes with a set of qualitative properties
        - Record the link properties {(node_id, node_id): link (neo4j object)}
        - Record the node properties {node_id: node (neo4j object)}
        - Let the crawl run along the edges until:
            - either the allowed number of steps to crawl is exhausted
            - there is no more nodes to use as a seed
    - DONE: change the weight calculation that will be using the link properties that were recorded
        - use the properties of the link and the node pair to calculate the weights for both
        matrices


 - NOPE: [FEATURE] [USABILITY] upon organism insertion and retrieval, use the 'orgnism' flag on the
        proteins and relationships to allow for simultaneous loading of several organisms.
        - Superseeded by the better way of doing it through multiple databases in a single neo4j
        instance

 - DONE: record the origin of the nodes and relationships:
    - Reactome
    - UNIPROT

 - DONE: define trust into the names of different databases and make use it as mask when pulling
    relationships

 - NOPE: [REFACTOR] remove the GraphDeclarator.py and re-point it directly into the cypher_drivers
    - It's already an abstract interface that can be easily re-implemented

 - NOPE: [REFACTOR] wrap the cypher_drivers into the db_io_routines class
    - Nope, it's already an abstract interface

 - DONE: [FUTUREPROOFING] [CODESMELL] get away from using `_properties` of the neo4j database
        objects.
        => Basically, now this uses a Node[`property_name`] convention

 - DONE: [PLANNED] implement the neo4j edge weight transfer into the Laplacian
    - DONE: trace the weights injection
    - DONE: define the weighting rules for neo4j
    - DONE: enable neo4j remote debugging on the remote lpdpc
    - DONE: change the neo4j password on remote lpdpc
    - DONE: add the meta-information for loading (eg organ, context, ...)
        - doable through a policy function injection

The other next step will be to register the context in the neo4j network in order to be able to
perform loads of networks conditioned on things such as the protein abundance in an organ or the
trust we have in the existence of a link.

    - neo4j database:
        - REQUIRE: add context data - basically determining the degree of confidence we want to
            have in the node. This has to be a property, because the annotations will be used as
            weights for edge matrix calculations and hence edges need to have them as well.

    - database parsing/insertion functions:
        - REQUIRE: add a parser to read the relevant information from the source files
        - REQUIRE: add an inserter to add the additional information from the parse files into
            the neo4j database
        - REQUIRE: an intermediate dict mapping refs to property lists that will be attached to
            the nodes or edges.

    - laplacian construction:
        - REQUIRE: a "strategy" for calculating weights from the data, that can reason on the in
            and out nodes and the edge. They can take in properties returned by a retrieval pass.
        - REQUIRE: the weighting strategy should be a function that can be plugged in by the end
            user, so for a form neo4j_node, neo4j_node, neo4j_edge > properties.
        - REQUIRE: the weighting strategy function should always return a positive float and be
            able to account for the missing data, even if it is raising an error as a response to
            it.
        - REQUIRE: the current weighting strategy will be encoded as a function using node types
            (or rather sources).


<DONE: CONFIGS sanity>

Current architecture:
    - `user_configs.py`, containing the defaults that can be overriden
    - `XXX.yaml` in  `$BIOFLOWHOME/configs` that allow them to be overriden
    - the override is performed in the `user_configs.py`
    - `main_configs` imports them and performs the calculation of the relative import paths
    - all other modules import `main_configs as configs` and use variables as `configs.var`

    - the injection is performed by an explicit if-else loop after having read in a sanitized way
    the yamls needed.

    - example_configs.yaml are deployed during the installation, that the user can modify.
    - the import will look for user_configs.yaml, that the user would have modified
    - the yaml file would contain a version of configs, se

    - DONE: how do we find where is $BIOFLOWHOME to read the configs from?
        - Look up the environment variable. If none found, proceed with default location
            (~/bioflow). ask user to register it explicitly upon installation.

    - DONE: define the copy to the user directory without overwrite

    - DONE: test that the edits did not break anything by defining a new $BIOFLOWHOME and pulling
        all online dbs again

    - DONE: move user_configs.py to the .yaml


A possible approach is to use the `cfg_load` library and the recommended good application practices:
    - define the configurations dictionary ()
    - load the defaults (stored in the configs file within the library)
    - load the $BIOFLOWHOME/configs/<> - potential overrides for the defaults
    - update the configurations with user's overrides
    - update the configurations with the environment variables (or alternatively read for the
        command line and inject into the environment)

    - PROBLEM: we use typing/naming assist from IDE
        We can override this issue by still assigning all the the values retrieved from the
    config files to the variable names defined in the main_configs file.
        In this way, our default variable definitions are the "default", whereas the user's yaml
    file supplies potential overrides

    - PROBLEM: storing the build/background parameters for the Interactome_Analysis and the
        Annotome_Analysis classes

Itermediate problem: there is a loading problem for the `BioKnowledgeInterface` due to `InitSet`
used on the construction (~6721 nodes) is significantly bigger than the `InitSet` used in order
to generate the version that is being put into storage.
    - `InitSet` is loaded as `InteractomeInterface.all_uniprots_neo4j_id_list`. In case
    `reduced_set_uniprot_node_ids` is defined in the parameters, the source set from the
    `InteractomeInterface` is trimmed to only nodes that are present in the limiter list.
    - `InteractomeInterface` is the one currently stored in a way that can be retrieved by a
    `.fast_load()` and is not changed since the last  since the last load. The `.fast_load()` call
    is performed

    - It looks like the problem is in the fact that the rebuild uses a background limiter (all
    detectable genes), whereas the fast load doesn't.

    => Temporary fix:
        - DONE: define a paramset with background in the `analysis_pipeline_example`.

    => In the end, it is a problem of organization of the parameters and of the context.
        - TODO: Set a user flag to know if we are currently using a background.
        - TODO: Set a checkers on the background loads to make sure we
        - TODO: Version the builds of the MolecularInterface and BioKnowledge interface
        - TODO: Provide a fast fail in case if the environment parameters differ between the
            build and the fastload. Environment parameters are in the `user_configs.py` file.
        - TODO: this is all wrapped in the environment variables
        - TODO: for each run, this is saved as .env text file render.


 - DONE: [SANITY] Configs management:
    - DONE: move the organism to the '~/bioflow'
    - DONE: all the stings `+` need to be `os.path.join`.
    - DONE: active organism is now the only thing that is saved. It is stored in "shelve" file
    inside the ".internal" directory
    - NOIMP: fold in the sources for the databases into a single location, with a selector from
    "shelve" indicating which organims to load.
    - NOIMP: create a user interface command in order to set up the environment and a saving file
    that allows the configs to be saved between the users.
    - DONE: move the `online_dbs.ini`, `mouse.ini`, `yeast.ini` to the `~/bioflow/configs` and add
    `user_configs.ini` to it to replace `user_configs.py`.

 - DONE: [SANITY]: move the location from which the base folder is read for it to be computed
        (for relative bioflow home insertion) (basically the servers.ini override)

The next step will be to register user configurations in a more sane way. Basically, it can
either be a persistent dump that is loaded every time the user is spooling up the program or an
.ini file in additions to the ones already existing

    - Peristent dump:
        - PLUS: Removes the need to perform reading in and out of .ini files
        - PLUS: Guarantees that the parameters will always be well formed
        - MINUS: is non-trivial to modify for the users

    - .ini files: => selected, but in the .yaml incarnation
        - PLUS: Works in a way that is familiar to most people
        - PLUS: Allows
        - MINUS: in case configurations are not properly defined, all crashes

    - Both:
        - NOPE: a command line process to define the variables
        - NOPE: a command line to show all the active flags
        - DONE: the configs folder to be gutted of active configs and them moved to the
            ~/bioflow directory
        - DONE: a refactor to show transparently the override between the default parameters
            and the user-supplied parameters
        - DONE: a refactor to remove the conflicting definitions (such as deployment vs test
            server parameters)


 - DONE: [USABILITY] add a proper tabulation and limit float length in the final results print-out
    (tabulate: https://pypi.python.org/pypi/tabulate)

 - DONE [USABILITY] add limiters on the p_value that is printed out elements

 - DONE: [USABILITY] change colors of significant elements to red; all others to black (with alpha)
        - This modification is to be peformed in the `samples scatter and hist` function in the
        `interactome_analysis` module

 - DONE: [FEATURE] [REFACTOR]:
    - add the selection of the degree width window for stat. significance calculation.

 - DONE: [FEATURE]:
    - add p-value and pp-value to the GO annotation export

 - DONE: Currently, performing an output re-piping. The output destinations are piped around thanks
    to a NewOutput class in main_configs, which can be initialized with a local output directory
    (and will be initialized in the auto-analyze function for both the interactone and the knowledge
    analysis network)

 - DONE: [DEBUG] add the interactome_network_stats.png to the run folder

 - DONE: [USABILITY]: fold the p-values into the GO_GDF export in the same way we do it for the
        interactome

 - NOIMP: [USABILITY] Add an option for the user to add the location for the output in the
        auto-analyse

 - DONE: [USABILITY] align the rendering of the conditions in annotations analysis with the
        interactome analysis

 - DONE: [DEBUG]: align BioKnowledgeInterface analysis on the InteractomeAnalysis:
    - Take the background list into account
    - Take in account the analytics UP list in the hashing (once background is taken into account)
    We are dealing with a problem on the annotation analysis network not loading the proper
    background (proably due to the wrong computation of the laplacian). At this point we need to
    align the Annotation analysis on the molecular analysis.
        - DONE: run git blame on the Molecular network interface, copy new modifications
        - DONE: run git blame on molecular network analysis, copy the new modifications

 - DONE: [SHOW-STOPPER]: debug why GO terms load as the same term
    - Not an issue - just similar GO terms of different types (eg entities)


 - TODO: [SANITY]: Feed the location of the output folders for logs with the main parameters
    - DONE: create a function to generate paths from a root location
    - DONE: define new "TODO"s : (TRACING, OPTIMIZE and CURRENT)
    - DONE: Move the "info" log outputs to the parameters
    - NOIMP: Allow the user to provide the names for the locations where the information will be
        stored
    - DONE: trace the pipings of the output / log locations

 - DONE: [USABILITY] add a general error log into the info files

 - DONE: [USABILITY] save the final table as a tsv into the run directory

 - DONE: [USABILITY] format the run folders with the list send to the different methods

 - DONE: [USABILITY] add a catch-it-all for the logs

 - DONE: [OPTIMIZATION]: Profile the runtime in the main loop:
    - DONE: check for consistency of the sparse matrix types in the main execution loop
    - DONE: run a profiler to figure out the number of calls and time spend in each call. common
        profilers include `cProfile` (packaged in the base python) and `pycallgraph` (although no
        updates since 2016). Alternatively, `cProfile` can be piped into `gprof2dot` to generate
        a call graph
    - DONE: biggest time sinks are:
        - csr_matrix.__binopt (likely binaries for csr matrices) (22056 calls, 81404 ms)
        - lil_matrix.__sub__ (7350/79 968)
        - lil_matrix.tocsr (7371/76 837)
        - sparse_abs (7350/58 700)
        - lil_matrix._sub_sparse (3675/48 192)
        - csr_matrix.dot/__mul__/_mul_sparse_matrix (~7350 / ~48 000)
        - triu (3675/47 592)
        - csr_matrix.tocsc (7353 / 47 253)
        - csc_matrix.__init__ (121332 / 47 253) (probably in-place multiplication is better)
    - DONE: first correction:
        - baseline: edge_current_iteration: 3675 calls, 362 662 ms, 40238 own time.
        - DONE: uniformize the matrix types towards csc
            - eliminated lil_matrix: performance dropped, lil_matrix still there)
                (3675 calls, 366 321ms, 41 525 own time)
            - elimitated a debug branch in get_current_through_nodes => No change
                (3675 / 367 066 / 40 238)
            - corrected all spmat.diag/triu calls to return csc matrices + all to csc => worse
                (3675 / 408 557 / 40 657)
            - tracked and imposed formats to all matrix calls inside the fast loop => 50% faster
                (3675 / 262 627 / 40 420)
                => csr_matrix still gets initialized a lot and a coo_matrix is somewhere.
                lil_matrix is gone now though
            - repaced all mat.csc() conversions by tocsc() calls
                => was more or less already done
    - DONE: profiled line per line execution.
        - sparsity changes are the slowest part, but seem unavoidable
        - followed by triu
        - followed by additions/multiplications
        - followed by cholesky
    - DONE: delay the triu until after the current accumulator is filled
        - baseline: 261 983
        - after: 197 888 => huge improvement
    - NOPE: perform in-place multiplication => impossible (no in-place dot/add/subract
        versions)
    - TODO: clean-up:
        - DONE: remove the debug filters connectors
        - DONE: deal with the confusing logic of enabling the splu solver
            - problem => We run into the optimization of using a shared solver
        - DONE: test the splu and the non-shared solver branch
            - Non-sharing works, is slow AF
            - Splu switch works, but is slow AF

 - DONE: resolve the problem with "memoization" naming convention. in our cases it's remembering
        potential diffs. It appears that it also interacts with "fast load" in "build extended
        conduction system. Technically, by performing a memoization into a database, we could
        have a searchable DB of past runs, so that the comparison is more immediate. So far the
        usage is restricted to `InteractomeInterface.compute_current_and_potentials`, to enable a
        `fast load` behavior

 - DONE: [DEBUG] [SHOW-STOPPER]: connections between the nodes seem to have disappeared
    - DONE: check if this could have been related to memoization. Unlikely. The
        only place where the results of memoized were accessed was for voltages > it is.
    - DONE: run test on the glycogen set. No problem detected there

 - DONE: [SANITY] Logs:
    => DONE: Pull the logs and internal dumps into the ~/bioflow directory
    => IGNORE: Hide away the overly verbose info logs into debug logs.

 - PTCH: [SANITY] allow user to configure where to store intermediates and inputs/outputs
 - DONE: [SANITY] move configs somewhere saner: ~/bioflow/ directory seems to be a good start
 - PTCH: [CRICITAL] MATPLOTLIB DOES NOT WORK WITH CURRENT DOCKERFILE IF FIGURE IS CREATED =>
        figures are not created.
 - DONE: [CRITICAL] ascii in gdf export crashes (should be solved with Py3's utf8)

 - DONE: [DEBUG]/[SANITY]: MongoDB:
    - DONE: Create a mongoDB connection inside the fork for the pool
    - DONE: Move MongoDB interface from configs into a proper location and create DB type -
            agnostic bindings

 - DONE: [SHOW-STOPPER] Memory leak debugging:
    - DONE: apply muppy. Muppy did not detect any object bloating > most likely comes from matrix
            domain
    - DONE: apply psutil-based object tracing. The bloat appears around summation + signature
            change operations in the main loop > sparse matrix summation and type change seem to be the origin.
    - DONE: try to have consistent matrix classes and avoid implicit conversions. Did not help
            with memory, but accelerated the main loop by 10x. Further optimization of the main loop might be
            desirable
    - DONE: try to explicitely destroy objects with _del and calls of gc. Did not mitigate the
            problem. At this point, the memory leak seem to be localized to C code for the
            summation/differentiation between csc_matrices.
    - DONE: disable multithreading to see if there is any interference there. Did not help
    - DONE: extract he summation to create a minimal example: did not help
    - DONE: build a flowchart to see all the steps in matrices to try to extract a minimal
            example. Noticed that memoization was capturing a complete sparse matrix. That's where the bloat
            was happening.
    - DONE: correct the memoization to remember the currents only.

 - DONE: Threading seem to be failing as well.
        The additional threads execute the first sampling, but never commit. Given that they
        freeze out somewhere in the middle, the most likely hypothesis is that they run out of
        RAM and only one thread - that keeps a lock on it  - continues going forwards.
        In the end it is due to memory leak.

 - DONE: [DEBUG]: sampling pools seem to be sharing the random ID now and not be parallel. CPU
        usage however indicates spawned processes running properly

 - DONE: Random ID assignements to the the threads seem to be not working as well
    - DONE: rename pool to thread
    - DONE: add the ID to the treads
    - DONE: debug why objects all share the same ID across threads (random seed behavior? - Nope).
            The final reason was that thread ID was called in Interactome_instance initialization and not

 - DONE: [USABILITY] change ops/sec from a constant to the average of the last run (was already
        the case)

 - DONE: debug the issue where the `all_uniprots_id_list` interesection with `background` leads
        to error-prone behavior. Errors:
        a) injection of non-uniprot IDs into the connected uniprots
        b) change of the signature of the Interface instance by changing:
            - `analytics_uniprot_list` in `InteractomeInterface`
            - `analytics_up_list` in `BioKnowledgeInterface`
        The issue seem to be stemming from the following variables:
        in `InteractomeInterface`:
            - `self.all_uniprots_neo4j_id_list` (which is a pointer to `self.reached_uniprots_neo4j_id_list`)
            - `self.connected_uniprots`
            - `self.background`
            - `self.connected_uniprots` and `self.background` are directly modified from the `auto_analyze`
                routine and then
            - The operation above is cancelled by random_sample specifically
                Which is probably the source of our problems. now the issue is how to get rid of the
                problem with nodes that failed
            - the issue only emerges upon sparse sampling branch firing
            - `self.entry_point_uniprots_neo4j_ids` is used by `auto_analyze` to determine sampling depth
                and is set by the `set_uniprot_source()` method and is checked by the
                `get_interactome_interface()` method
        in `BioKnowledgeInterface`:
            - `self.InitSet` (which is `all_uniprots_neo4j_id_list` from the `InteractomeInstance` from
                which the conduction system is built)
            - `self.UPs_without_GO`
        The intermediate solution does not seem to be working that well for now: the sampling mechanism
        tends to pull as well the nodes that are not connected to the giant component in the neo4j graph.
        - Tentatively patched by making the pull from which the IDs are sampled stricter. Seems to
        work well

 - DONE: [SHOW-STOPPER]: ReactomeParser does not work anymore likely a node.
        The issue was with a automated renaming during a refactoring to extract some additional
        data

 - DONE: [FEATURE]: (done by defining a function that can be plugged to process any tags in neo4j)
    - In Reactome, parse the "Evidence" and "Source" tags in order to refine the laplacian weighting


Bigger refactors
****************

 - TODO: [OPTIMIZATION]:
    - bulk-group the insertions and cross-linkings for the Reactomse

 - TODO: [USABILITY] pull inlined updates printing from evoGANs project.
    => Currently the percentages are managed by log.info(calls)
    => Providing an in-line update would require a print(<log message>, end='\r')
    => Change log management so that the info gets logged into a file without rising to the
        surface and couple all the log.info with a "print"

 - TODO: [FEATURE] [DECIDE IF IMPLEMENT] add flow calculation for real samples saving to mongodb, +
buffering (if unchanged InteractomeInstance and other secondary formatting, just retrieve the
flow from the database)

 - TODO: [REFACTOR] refactor the entire edge typing upon insertion, retrieval upon construction of
laplacian/adjacency matrix and setting of Laplacian weights

 - TODO: [FEATURE] Change the informativity between nodes connection from the max of their
informativities
to the difference of their informativities. In this way, the total path is equivalent to the
quantity of the information stored

 - TODO: [FEATURE] provide the interface for overlaying the molecular maps to check for
signatures/compare samples

 - TODO: [FEATURE] change the informativity computation so that the path between nodes through the
network is log of probability of being connected through that suite of GO terms.

 - TODO: [FEATURE] flow to a targeted set

 - TODO: [FEATURE] we ighted targets flow

 - TODO: [FEATURE] modification of the Laplacian weights by the end user.

 - TODO: [FEATURE] import credence from the interaction databases

 - TODO: [FEATURE] add direct interactions for TFs in yeasts by combining
https://www.nature.com/articles/ng2012 and https://pubmed.ncbi.nlm.nih.gov/29036684/ for binary
interaction filters

 - TODO: [DEV TOOL] performance evaluation run: compute compops for sampling a large pool of
genes in yeast

 - TOOD: [REFACTOR] set the model for the database usage for the samples by performing
insertions/conversions inside the `sample_storage` databases wrappers (for mongodb: form the dict
of query/payload)

 - TODO: [NOT THIS ITERATION] neo4j and mongodb versionning based on the currently
active organism. For instance, all neo4j nodes and edges need to be marked with the organism tag.

 - TODO: [NOT THIS ITERATION] add a mail signalling to indicate the termination or crash of the
execution



Bulk Backlog:
-------------

Language of network alignment/explanation of net1 by net2: allows
   to compare GO annotation with interactome, cell type specificity
   analysis or organ context.
    => Solved by the search of average resistance between connected nodes in the graph1 on graph2
   . Alternatively, average flow intensity difference between a random sample of common noces on
   graphs 1 and 2.

Problems uncovered with user while testing the Docker integration:
    - explain how directories on the OS are mapped to the directories on the Docker
    - Suggest that in case you are using Mac/OSX, you need to manually increase memory allocated to Docker to at least 16 GBs:
        - 2 GB for each database
        - + ~7GB for each processor used to perform random sampling

Potential problems with pip installation:
    - the configs/dump files modified by the user
        => Move them to the ~/bioflow/ directory

Travis tester:
    - unittest
    - docker from Githubs
    - pip install
    - docker-compose
    - databases downloads

=> Build an export of the sampling current weights to figure out which nodes are offending.
=> Re-compute eigenvalues of the laplacian and add values to the network weighting nodes

=> Modify the insertion/retrieval pattern into the mongoDB to separate foreground from background runs
    - Start with the foreground run
    - continue with the backgrounds until satisfied
    - always possible to resume the sampling afterwards

Add a mention for what were the parameters of the launch of the analysis - what was build and where the data was loaded from?

We are using Interactome Interface for 5 independent reasons:
    - build the laplacian matrix
    - store the laplacian matrix
    - perform sampling on the laplacian matrix
    - calculate the stats on teh sampling matrix => This is actually done in interactome analysis
    - (OK) export the rendering to the gdf => this actually is done by a separate object.
    - it might be a good idea to refactor it.

We can already factor out the two methods responsible for a laplacian matrix building.

(OK) Correct the HINT downloading and renaming

Switch to matrix instead of dict for a current/tension storing in a dense fashion

(OK) Implement output redirects - main_config Outputs patching does not seem to work - we need
to create the object anew in case of need.

# correct the overestimation of flow gain for low-edge nodes in the network.

ADD to the documentation:
    - management of multicast to accession numbers and gene names.
        - random returns
        - cross-linked with 'is_likely_similar' links, that are imported to Laplacian with

Logging and CLI wrappers:
    - redirect logging to stderr
    - add version flag (version + commit #)
    - add autocomplete
    - dump gdf to stdout?
    - check option prompts
    - provide an interface to inform of the program completion (?)
    - add a spinner for slow processes

Documentation:
    - installation:
        - docker
        - Ubuntu
    - usage:
        - core command lines
        - core Python
        - post-processing for analysis
    - neo4j manual usage
        - Cypher
        - access to non-local neo4j instance
        - useful commands
        - what to do if the commands are slow (optimized for use case of the Bioflow, not necessarily best)


Bulk Backlong Done:
-------------------

YEAST:
    DISABLE TRRUST/human/mouse-specific imports: DONE

Next steps, in order:
    - (DONE) dump of indexed nodes Legacy Ids and a method to compare them (in the "DB inspection" realm)
    - (DONE) Delete dead branches, break dependency on bulbs
        - think if we could do testing for a neo4j build
    - (DONE) build a new docker image

There still seem to be a problem of regular convergence to the same paths in the network. Potential sources:
    - borked topology
    - current intensity between the interconnected nodes (potentially resolved)
    - tight clusters due to cross-linking that disperse the network
=> Solved through statistics computation refactoring


Functional:
-----------
-   (OK )Enforce p-value limitation to what is achievable with the background sampling size.

-   (OK ) Re-inject the p_values following the comparison back into the output and export them as part of
    gdf schema

-   (OK - Won't fix) Normalize the Laplacian (Joel Bader)

-   (OK - Won't fix) Increase the fudge factor to ~ 10% of error (Joel Bader)

-   (OK) Integrate the amplification level into the analysis - relative amplitude of perturbation

    -   Would require to fold the hits list into the interface object

    -   Would require to always use the hits list as

    -   Would require an explicit injection of voltage (rather a normalization)

-   For dense sampling, switch to a matrix instead of a dictionary to store current/tension for each
    pair of proteins.

-   Enable a calculation with explicit sinks and sources groups

    -   Would require a split in hits into "sink" and "source" groups

    -   Would need a revised null model (random sinks/sources? if single sink, random sources?)

-   (OK) Add signal/noise ration - the flow we are getting in a given node compared to what we would have
    expected in a random node.

        -   => More or less already done by the p_value; excpet p_value also accounts for the node of average degree X

-   Wrap neo4j and laplacian files loading/offloading into top-level commands

    -   check if this is compatible with long-term neo4j architecture => have to; multi-org systems

-   Flow system comparison

    -   change in the Laplacian

    -   comparison of specific flow patterns

    -   separation between instance, actual flow and flow comparator.

Structural:
-----------
-   (OK) Create a separate structure for performing the statistical analysis, that is independent from the
    wrapper

-   (OK) Enforce the single source of the Interface Objects for sampling to simplify
    consistency enforcing

-   Return the control of the reach-limiter from the Knowledge interface to the current routines (Why?)

-   (DONE) Insert the effective sample into the storage DB to avoid re-running it in case more background is
    needed or the background needs to be interrupted. (Implemented otherwise)

-   (DONE) In the random reference generator, reformat so that different types of queries are matched
    with the same types of queries. This is required in order to provide support for statistics on
    multiple query types calls.

    -   There are two levels of tagging: conduction system (eg background, ....) and query system
        (hits, circulation, etc...)

    -   There is also a need for specific querying for the conduction system to avoid re-generating
        them

    -   As well as "hit" analysis systems, so that they can be retrieved instead of needing to be
        re-generated.

-   (Won't fix)Move node annotation loading/offloading to an elasticsearch instance, always mapping to uniprots

-  (OK) The alteration of the chain of statistics calculation is somewhat hardcore

Integration:
------------
- (DONE) Update to CYPHER and a more recent neo4j instance accessed through bolt
    - possibility of using a periodic Cypher update (each 500-1000) instead of atomic?

- (Rejected) Consider removing neo4j altogether => Nope, it's a good persistence solution


Cosmetic:
---------

-   Properly indent multi-line :param <parameter type> <parameter name>: descriptors

-   (OK) Integrate compops/second estimation to the sources.ini

-   Perform profiling by creating a dedicated set of loggers that would log an "execution time"
    flag set

-   Add the forwarding of the thread number to the progress report

-   (CLI) Reformat progress report as a progress bar.

-   (DONE) The explanation of the model and its transmission into the published data needs to be easy enough
    to explain why the hits are justified and why they aren't.


Feature wishlist:
-----------------

-  Add protein abundance level for instantiation of the network

-  Add a coarseness feature on the interactome analysis affecting
   sampling behavior, so that precision is sacrificed in favor of
   computation speed.

-  Build an "inspector routine" that would allow us to see the nodes that would allow us to route
   the most information => we need to recompute the most central nodes in Interactome, because we
   still observe a heavy skew in the nodes with a high degree.

-  We always need to first build Interactome Interface before BioKnowledge Interface and in the
   end we need to have both of them build before we can run automated analysis. A nice fix would be
   to raise flags when they are loaded, instead of relying on the loading behaviors. (Deprecated)

-  (WillNotFix) Event sourcing pattern for the graph assembly and modification from the base databases. (Far Fetched)

-  The execution entry points have to be the five canonical queries.

-  Write the flow groups so that it is possible to calculate the information circulation
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
   into a single database (Long-term applicability enhancement).

-  Clustering algorithms going beyond the spectral clustering,

    - Not needing a pre-defined number of clusters

    - Able to assign the same node to several clusters

    - Maybe iterative DBSCAN or agglomerative clustering with removal of detected clusters until
      we hit some threshold on the number of number of nodes - average circulation in cluster curve
      obtained from random nodes sampling

    - We can deduce that graph from the clustering of sets of random nodes v.s.

-  Graph exploration module:

    - Strongest eigenvectors / highest circulation in a random set of nodes

    - Randomized clustering

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


-  (DONE) separate the envelopes for the GO and Reactome graphs retrieval from
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

-  Inline the background for the InteractomeInstance into the __init__

-  Inline the undumping and dependent variables calculation into the __init__ of
   InteractomeInstance and BioKnowledgeInterface

-  change to element import directories from which too many
   functions/objects are imported (import numpy as np)

-  Make methods running large systems of procedures to being
   dictionary-driven (Why?)

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
   unittests) will switch it to a mock, not expecting the database to answer

-  (DONE) Bulk-insertion into the neo4j. => Requires taking over the bulbs engine

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

-  (DONE) build a condas-compatible package that would be installable
   cross-plateform and would contain pre-compiled binaries for
   C-extenstions. (Failed - we are better off with a Docker given complexity of the stack)

-  In case we are calling time-consuming parsers from several location,
   we might want to insert "singleton" module into the block, that
   performs all the parsings only once per program run.

-  We need to dynamically update the values of main_config whenever the location whenever
   configs from configsfiles are modified, so their modification do not require restaring the
   program. Alternatively we say that the configs need to be modified before the rest of the
   program can be executed.

-  Transform all the matrices so that the first one is packed line-based and the second one
   column-based. This will allow the optimization for register pre-pulling in the processor

-  Docker container cp command to accelerate the database rebuild process?

Neo4j and Laplacian construction:
---------------------------------

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


New Databases:
--------------

- Protein abundance

- (DONE) Transcription/translation regulation

- Post-translational modifications

- Isoforms

- Interactions that store annotation datasets, such as IntAct


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

- add pictures of what netqork analysis looks like

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
   input matrices we would like to inspect and what output matrices we
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
   many null eigenvalues as there are connex segments in the graph (Error is acceptable, LU is slower)

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

-  TODO: Add the 95% confidence interval for a given precision rate for the depth of sampling.
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

-  Database sources cited in the differential network biology paper by {Idecker 2012}:
    - BioGRID
    - HPRD (Human Reference Protein-protein interaction Dataset - human only)
    - IntAct (good idea to integrate given the quality and extensivity of data)
 