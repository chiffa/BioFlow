"""
Saner io and filesystem manipulation compared to python defaults
"""
import os
from shutil import rmtree


def mkdir_recursive(path):
    """
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param path:
    :return:
    """
    print 'trying to read: ', path
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    print 'subpath: %s' % directory_name
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(path):
        print 'path %s does not exist yet, looks like file: %s' % (path,
                                                                   '.' in path.split('/')[-1][-5:]),
        if '.' not in path.split('/')[-1][-5:]:  # should be able to suppress specific file creation
            os.mkdir(path)
            print '; created'
        else:
            print '; creation skipped'


def wipe_dir(path):
    """
    wipes the indicated directory
    :param path:
    :return: True on success
    """
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    print directory_name
    if not os.path.isdir(path):
        print Exception('failed to delete %s: not a directory')
        return False
    if not os.path.exists(directory_name):
        return True  # Nothing to do: destruction already done
    else:
        rmtree(path)
        return True
