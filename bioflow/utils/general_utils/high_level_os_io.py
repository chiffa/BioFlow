"""
Saner io and filesystem manipulation compared to python defaults
"""
import os
from shutil import rmtree, copyfile
from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def mkdir_recursive(path):  # pragma: no cover
    """
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param path:
    :return:
    """
    log.debug(
        'trying to create recursively path containing: %s', path)
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    log.debug('subpath: %s', directory_name)
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(path):
        log.debug('path %s does not exist yet, looks like file: %s',
                  path, '.' in path.split('/')[-1][-5:])
        if '.' not in path.split('/')[-1][-5:]:
            # should be able to suppress specific file creation
            os.mkdir(path)
            log.debug('; created')
        else:
            log.debug('; creation skipped')


def wipe_dir(path):  # pragma: no cover
    """
    wipes the indicated directory
    :param path:
    :return: True on success
    """
    path = os.path.abspath(path)
    log.debug('entered ')

    if not os.path.exists(path):
        log.debug('path does not exist')
        return True  # Nothing to do: destruction already done

    if os.path.isdir(path):
        directory_name = path
    else:
        directory_name = os.path.dirname(path)

    log.debug('going to wipe %s for path %s', directory_name, path)
    if not os.path.isdir(directory_name):
        log.exception(
            'failed to delete %s: for path %s, not a dir', directory_name, path)
        return False

    for sub_path in os.listdir(directory_name):
        if os.path.isdir(sub_path):
            log.exception(
                'failed to delete %s: for path %s, anti rm -rf flag', directory_name, path)
            return False

    log.debug('performing a rmtree')
    rmtree(directory_name)
    return True


def copy_if_doesnt_exist(source, destination):  # pragma: no cover
    """
    Copies a file if it does not exist

    :param source:
    :param destination:
    :return:
    """
    if not os.path.isfile(destination):
        log.info('Did not detect file %s at the destination, copying from source %s' % (source,
                                                                                        destination))
        copyfile(source, destination)