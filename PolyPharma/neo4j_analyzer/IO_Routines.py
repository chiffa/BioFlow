__author__ = 'ank'

import pickle


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
    pickle.dump(object_to_dump, DF)
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
    return pickle.load(DF)