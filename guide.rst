Assembling the files required for the database creation:
========================================================

In order to build the database, the program is going to look for the following files specified
in the following locations within the PolyPharma/configs/sources.ini::

    * OBO 1.2 file of GO terms and relations, available at: http://www.geneontology.org/GO.downloads.ontology.shtml
    * will look for at the path in [GO] - "location"

    * UNIUPROT-SWISSPROT .txt text database file available at: http://www.uniprot.org/downloads
    * will look for the files at the path [UNIPROT] - "location"
    * will load the information for the organism with specified NCBI taxonomy identifier from "tax_id"

    * Reactome.org "Events in the BioPax level 3" file, available at: http://www.reactome.org/download/index.html
    * will look for the files at [REACTOME] - "location"
    * will only load the file specified by the "load" parameter

    * HiNT binary interaction files for the organisms of interest, availble at: http://hint.yulab.org/batch.html
    * will look for the files at the path [HINT] - "location"
    * will load the information for the organism with specified NCBI taxonomy identifier from "load" parameter

    * BioGRID ALL_ORGANISMS file in the tab2 format, available at http://thebiogrid.org/download.php
    * will look for for the files at the path [BIOGIRD] - "location"
    * will load only the file specified in the "load" parameter

    * Gene-chromosome mapping files from the Uniprot documentation: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ (needed only for working on aneuploidy)
    * Configured in the [CHROMOSOMES] section
    * currently not ready for public use

    * Organism-specific protein aboundance files, available at: http://pax-db.org/#!downloads
    * Configured in the [ABOUNDACES] section
    * currently not ready for public use


It is possible to specify the file locations and indentifiers manually, and then download and install them
to the specified locations manually.

However the following command should be able to do it for you for three commonly used organism (human, mouse, saccharomyces cerevisae),
provided you follow the instructions properly::

    > python CLUI.py initialize --path myfolder --neo4jserver http://localhost:7474 --mongoserver mongodb://localhost:27017/

    > python CLUI.py downloaddbs

    > python CLUI.py setorgconfs --organism [mouse, human, yeast]

Typical sources.ini configfile:
===============================

Here is what a typical configfile would look like::

    [REACTOME]
    location = myfolder/External_DBs_Store/Reactome
    load = Mus musculus.owl

    [UNIPROT]
    location = myfolder/External_DBs_Store/Uniprot/uniprot_sprot.dat
    tax_ids = 10090,

    [HINT]
    location = myfolder/E/External_DBs_Store/HiNT
    load = MouseBinaryHQ.txt

    [GO]
    location = myfolder/E/External_DBs_Store/GO/go.obo

    [BIOGIRD]
    location = myfolder/E/External_DBs_Store/BioGIRD
    load = Mus_musculus.tsv

    [CHROMOSOMES]
    location = myfolder/E/External_DBs_Store/Chr_mappings
    load = mouse
    namepattern = mouse

    [ABOUNDANCES]
    location = /home/ank/Documents/External_DBs_Store/Protein_aboundances
    load = 10090

The data relative to the following parameters::

    [REACTOME]

    [UNIPROT]

    [HINT]

    [GO]

    [BIOGIRD]

is critical for any application and must be properly configured and is critical for any application
of the method.

On the other hand the following parameters are here for legacy application reasons and are not currently
documented::

    [CHROMOSOMES]

    [ABOUNDANCES]

the "load" parameter in the "[UNIPROT]" folder requires a trailing comma and can take in multiple arguments
separaged by a comma and a space, in case UNIPROT identifiers of proteins from several organisms are desired
(for instance when host-disease proteome interactions are investigated)

The configuration files might be declared and switchedmanually (only the "source.ini" one will be parsed,
folders such as "sources_organism.ini" will be ignored and can be renamed to "source.ini" quite easily)

It is possible for the users to generate source.ini file for three organisms with the following command::

    python CLUI.py setorgconfs --organism [mouse, human, yeast]

This allows to switch rapidly between different investigated organism.

Please don't forget to switch or purge neo4j databases between organisms, because each organism needs it's own neo4j instance.

Database manipulation:
======================


In case you are going to work with organisms with large proteomes (mouse, human), neo4j might run out of memory and prompt to be restarted with
a larger allocation of RAM. Please, follow the instructions, then empty the database and start the import again
