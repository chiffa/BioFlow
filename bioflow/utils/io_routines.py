"""
Defines a couple of useful method to perform IO on dumping files
"""
from pickle import load, dump
from csv import reader
from bioflow.configs.main_configs import Dumps
from time import time
import subprocess
import numpy as np


def _is_int(_obj):
    """
    Checks if an object is an int with a try-except loop

    :param _obj:
    :return:
    """
    try:
        int(_obj)
    except TypeError or ValueError as e:
        return False
    else:
        return True


def _get_git_revision_hash():
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'])


def _get_git_revision_short_hash():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])


def write_to_csv(filename, array):
    """
    Dumps a numpy array to csv

    :param filename: location of dumping
    :type filename: str
    :param array: array to dump
    :type array: numpy.array
    """
    np.savetxt(filename, array)
    # dump_file = open(filename, 'wt')
    # dump_file.write(array)
    # dump_file.close()


def dump_object(dump_filename, object_to_dump):
    """
    Shortcut for pickling & dumping behavior

    :param dump_filename: filename where the object will be dumped
    :type dump_filename: str

    :param object_to_dump: object to be pickled and dumped
    :type object_to_dump: pickable object
    """
    dump_file = open(dump_filename, 'wb')
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
    dump_file = open(dump_filename, 'rb')
    # print(dump_filename)
    return load(dump_file)


def deprecated_get_bulbs_ids_set(location):
    """
    Retrieves bulbs ids for the elements for the analyzed group

    :param location: where the bulbs ids are loaded

    """
    bulbs_ids = []

    with open(location, 'rt') as src:
        csv_reader = reader(src)
        for row in csv_reader:
           bulbs_ids.append(row)

    print('debug: %s' % bulbs_ids)

    if len(bulbs_ids[0]) == 1:
        bulbs_ids = [int(ret[0]) for ret in bulbs_ids]

    else:
        bulbs_ids = [(int(ret), float(ret_w)) for ret, ret_w in bulbs_ids]

    return bulbs_ids


def get_source_bulbs_ids():  # TRACING: signature change from 1 to 2 returns
    """ retrieves bulbs ids for the elements for the analyzed group """
    return undump_object(Dumps.analysis_set_bulbs_ids)


def get_background_bulbs_ids():
    """ retrieves bulbs ids for the elements in the background group """
    return undump_object(Dumps.background_set_bulbs_ids)


def memoize(f):
    """"
    The standard memoization wrapper for a function
    """
    memdict = {}

    def internal_function(*args, **kwargs):
        if args in list(memdict.keys()):
            return memdict[args]
        else:
            result = f(*args, **kwargs)
            memdict[args] = result
            return result
    return internal_function


def time_exection(f):
    """
    The standard timing wrapper for a function
    """
    def int_function(*args, **kwargs):
        now = time()
        result = f(*args, **kwargs)
        print(time() - now)
        return result

    return int_function


if __name__ == "__main":
    # REFACTOR: [better environment]: actually inject the revision hashes into the execution
    #  pipeline in order to allow for reproductible execution
    print(_get_git_revision_hash())
    print(_get_git_revision_short_hash())
