__author__ = 'ank'

import pickle
from pickle import load, dump

def write_to_csv(filename, array):
    """
    Dumps a numpy array to csv

    :param filename: location of dumping
    :type filename: str
    :param array: array to dump
    :type array: numpy.array
    """
    DF = file(filename, 'w')
    DF.write(array)
    DF.close()


def dump_object(dump_filename, object_to_dump):
    """
    Shorcut for pickling & dumping behavior

    :param dump_filename: filename where the object will be dumped
    :type dump_filename: str

    :param object_to_dump: object to be pickled and dumped
    :type object_to_dump: pickable object
    """
    DF = file(dump_filename, 'w')
    dump(object_to_dump, DF)
    DF.close()


def undump_object(dump_filename):
    """
    Undumps a pickled object

    :param dump_filename: filename from which undump
    :type dump_filename: str
    :return: the undumped object
    :rtype: object
    """
    DF = file(dump_filename, 'r')
    return load(DF)

if __name__ == "__main":
    pass

# def dump_object_to_mongoDB(object, object_properties_dict):
#     """
#     A general interface to dump the objects to the mongo database
#
#     :param object: object we are willing to dump
#     :param object_properties_dict: string dict of properties of the object we are willing to dump
#     """
#     object_properties_dict.update({'Payload':pickle.dumps(object)})
#     tmp_coll.insert(object_properties_dict)
#
#
# def recover_and_undump_from_MongoDB(property_index_dict, upper_bound = 30):
#     """
#     A general interface to recover the objects from the mongo database based on a property list
#
#     :param upper_bound: upper bound on how many payload-containing packages the request could return
#     :param property_index_dict: string dict mapping the property keys to the property values
#     :return: list of object containers responding to the property
#     """
#     if tmp_coll.find(property_index_dict).count()>upper_bound:
#         return Exception('Number of obejcts covered by this cursor exceeds the upper bound specified in the paramenters')
#     postlist = []
#     for post in tmp_coll.find(property_index_dict):
#         postlist.append(post)
#
#     return postlist
