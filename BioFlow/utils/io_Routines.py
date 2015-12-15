"""
Defines a couple of useful method to perform IO on dumping files
"""
from pickle import load, dump


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

if __name__ == "__main":
    pass
