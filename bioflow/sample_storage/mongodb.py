from pymongo import MongoClient
from bioflow.configs.main_configs import mongo_db_url, pymongo_prefix, pymongo_suffix


# Original module ported from main_configs
client = MongoClient(mongo_db_url)
db = client.BioFlow_database
annotome_rand_samp = db[pymongo_prefix + "UP_r_samples" + pymongo_suffix]
interactome_rand_samp_db = db[pymongo_prefix + "Interactome_samples" + pymongo_suffix]


# This refactoring is done to avoid problems with forks in python multithreading. Given writes
# and reads are rare in our architecture, the performance penalty should be reasonable.
# TODO: the problem is that the client is still not getting closed, which might lead to problems
#  with the threads pooling. => with wrapper?

# with MongoClient(mongo_db_url) as client:
#     loc_annotome_rand_samp(client).drop()

def loc_annotome_rand_samp():
    client = MongoClient(mongo_db_url)
    return client.BioFlow_database[pymongo_prefix + "UP_r_samples" + pymongo_suffix]


def loc_interactome_rand_samp():
    client = MongoClient(mongo_db_url)
    return client.BioFlow_database[pymongo_prefix + "Interactome_samples" + pymongo_suffix]


def drop_all_annotome_rand_samp():
    loc_annotome_rand_samp().drop()


def drop_all_interactome_rand_samp():
    loc_interactome_rand_samp().drop()


def insert_annotome_rand_samp(payload_dict):
    loc_annotome_rand_samp().insert(payload_dict)  # => inset_one()


def insert_interactome_rand_samp(payload_dict):
    loc_interactome_rand_samp().insert(payload_dict)  # => inset_one()


def find_annotome_rand_samp(filter_dict):
    return loc_annotome_rand_samp().find(filter_dict)
    # => iterate through the cursor to close the connection


def find_interactome_rand_samp(filter_dict):
    return loc_interactome_rand_samp().find(filter_dict)
    # => iterate through the cursor to close the connection


def count_annotome_rand_samp(filter_dict):
    return loc_annotome_rand_samp().count_documents(filter_dict)


def count_interactome_rand_samp(filter_dict):
    return loc_interactome_rand_samp().count_documents(filter_dict)

