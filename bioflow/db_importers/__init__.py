"""
This module manages the import of database connections into the Neo4j instance. In case the
weights are presented as ID_1, ID_2, weight, the insertion can be done by adapting the "insert
into the database" method from almost any database insertion.

Pasring and insertion of the files in the BioPax format is more complex - please refer to the
Reactome parsing and inserting
"""