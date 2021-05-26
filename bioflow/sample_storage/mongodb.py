"""
thinly wrapped MongoDB backend for the sample storage
"""
from pymongo import MongoClient
from bioflow.configs.main_configs import mongo_db_url, pymongo_prefix, pymongo_suffix


# Original module ported from main_configs
client = MongoClient(mongo_db_url)
db = client.BioFlow_database
# REFACTOR: [Better database]: change mongoDB names to something more intuitive

annotome_rand_samp = db[pymongo_prefix + "UP_r_samples" + pymongo_suffix]
interactome_rand_samp_db = db[pymongo_prefix + "Interactome_samples" + pymongo_suffix]


# This refactoring is done to avoid problems with forks in python multithreading. Given writes
# and reads are relatively rare in our architecture, the performance penalty should be reasonable.

def loc_annotome_rand_samp():
    """loads a session for database connection for annotome samples"""
    client = MongoClient(mongo_db_url)
    return client.BioFlow_database[pymongo_prefix + "UP_r_samples" + pymongo_suffix]


def loc_interactome_rand_samp():
    """loads a session for database connection for interactome samples"""
    client = MongoClient(mongo_db_url)
    return client.BioFlow_database[pymongo_prefix + "Interactome_samples" + pymongo_suffix]


def drop_all_annotome_rand_samp():
    """ drops all annotome samples"""
    loc_annotome_rand_samp().drop()


def drop_all_interactome_rand_samp():
    "drops all interactome samples"
    loc_interactome_rand_samp().drop()


def insert_annotome_rand_samp(payload_dict):
    """
    Adds a sample from annotome run

    :param payload_dict:  sample contents
    :return:
    """
    loc_annotome_rand_samp().insert_one(payload_dict)


def insert_interactome_rand_samp(payload_dict):
    """
    Adds a sample from the interactome run

    :param payload_dict: sample contents
    :return:
    """
    loc_interactome_rand_samp().insert_one(payload_dict)


def find_annotome_rand_samp(filter_dict):
    """
    Finds a sample in the annotome database

    :param filter_dict: arguments dict according to which perform the search
    :return:
    """
    return loc_annotome_rand_samp().find(filter_dict)
    # => iterate through the cursor to close the connection


def find_interactome_rand_samp(filter_dict):
    """
    Finds a sample in the interactome database

    :param filter_dict: arguments dict according to which perform the search
    :return:
    """
    return loc_interactome_rand_samp().find(filter_dict)
    # => iterate through the cursor to close the connection


def count_annotome_rand_samp(filter_dict):
    """
    Number of samples in the annotome that satisfy the filtering conditions

    :param filter_dict: arguments dict according to which perform the filtering
    :return:
    """
    return loc_annotome_rand_samp().count_documents(filter_dict)


def count_interactome_rand_samp(filter_dict):
    """
    Number of samples in the interactome that satisfy the filtering conditions

    :param filter_dict: arguments dict according to which perform the filtering
    :return:
    """
    return loc_interactome_rand_samp().count_documents(filter_dict)

