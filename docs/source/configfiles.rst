Configuration Files Description:
================================

Typical sources.ini configuration file content::

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
[SIDER]
location = myfolder/E/External_DBs_Store/SIDER2/meddra_adverse_effects.tsv
[ABOUNDANCES]
location = myfolder/E/External_DBs_Store/Protein_aboundances
load = 10090
[CHROMOSOMES]
location = myfolder/E/External_DBs_Store/Chr_mappings
load = mouse
namepattern = mouse

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

[SIDER]
[ABOUNDANCES]
[CHROMOSOMES]

the "load" parameter in the "[UNIPROT]" folder requires a trailing comma and can take in multiple arguments
separaged by a comma and a space, in case UNIPROT identifiers of proteins from several organisms are desired
(for instance when host-disease proteome interactions are investigated)

The configuration files might be declared and switchedmanually (only the "source.ini" one will be parsed,
folders such as "sources_organism.ini" will be ignored and can be renamed to "source.ini" quite easily)

It is possible for the