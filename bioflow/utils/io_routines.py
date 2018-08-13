"""
Defines a couple of useful method to perform IO on dumping files
"""
from pickle import load, dump
from csv import reader
from bioflow.main_configs import Dumps
from time import time
import subprocess


def get_git_revision_hash():
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'])


def get_git_revision_short_hash():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])


def write_to_csv(filename, array):
    """
    Dumps a numpy array to csv

    :param filename: location of dumping
    :type filename: str
    :param array: array to dump
    :type array: numpy.array
    """
    dump_file = file(filename, 'w')
    dump_file.write(array)
    dump_file.close()


def dump_object(dump_filename, object_to_dump):
    """
    Shortcut for pickling & dumping behavior

    :param dump_filename: filename where the object will be dumped
    :type dump_filename: str

    :param object_to_dump: object to be pickled and dumped
    :type object_to_dump: pickable object
    """
    dump_file = file(dump_filename, 'w')
    dump(object_to_dump, dump_file)
    dump_file.close()


def undump_object(dump_filename):
    """
    Undumps a pickled object

    :param dump_filename: filename from which undump
    :type dump_filename: str
    :return: the undumped object
    :rtype: object
    """
    dump_file = file(dump_filename, 'r')
    return load(dump_file)


def get_bulbs_ids_set(location):
    """
    Retrieves bulbs ids for the elements for the analyzed group

    :param location: where the bulbs ids are loaded

    """
    bulbs_ids = []
    with open(location) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            bulbs_ids = bulbs_ids + row
    bulbs_ids = [int(ret) for ret in bulbs_ids]
    return bulbs_ids


def get_source_bulbs_ids():
    """ retrieves bulbs ids for the elements for the analyzed group """
    return get_bulbs_ids_set(Dumps.analysis_set_bulbs_ids)


def get_background_bulbs_ids():
    """ retrieves bulbs ids for the elements in the background group """
    return get_bulbs_ids_set(Dumps.background_set_bulbs_ids)


def memoize(f):

    memdict = {}

    def internal_function(*args, **kwargs):
        if args in memdict.keys():
            return memdict[args]
        else:
            result = f(*args, **kwargs)
            memdict[args] = result
            return result
    return internal_function


def time_exection(f):

    def int_function(*args, **kwargs):
        now = time()
        result = f(*args, **kwargs)
        print time() - now
        return result

    return int_function


if __name__ == "__main":
    print get_git_revision_hash()
    print get_git_revision_short_hash()
