"""
Saner io and filesystem manipulation compared to python defaults
"""
import os
from shutil import rmtree
from BioFlow.utils.LogManager import logger


def mkdir_recursive(path):
    """
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param path:
    :return:
    """
    logger.debug(
        'trying to create recursively path containing: {0}'.format(path))
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    logger.debug('subpath: {0}'.format(directory_name))
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(path):
        logger.debug('path {0} does not exist yet, looks like file: {1}'.format(
            path, '.' in path.split('/')[-1][-5:]))
        if '.' not in path.split(
                '/')[-1][-5:]:  # should be able to suppress specific file creation
            os.mkdir(path)
            logger.debug('; created')
        else:
            logger.debug('; creation skipped')


def wipe_dir(path):
    """
    wipes the indicated directory
    :param path:
    :return: True on success
    """
    path = os.path.abspath(path)
    logger.debug('entered ')

    if not os.path.exists(path):
        logger.debug('path does not exist')
        return True  # Nothing to do: destruction already done

    if os.path.isdir(path):
        directory_name = path
    else:
        directory_name = os.path.dirname(path)

    logger.debug('going to wipe {0} for path {1}'.format(directory_name, path))
    if not os.path.isdir(directory_name):
        logger.exception(
            'failed to delete {0}: for path {1}, not a dir'.format(directory_name, path))
        return False

    for sub_path in os.listdir(directory_name):
        if os.path.isdir(sub_path):
            logger.exception(
                'failed to delete {0}: for path {1}, anti rm -rf flag'.format(directory_name, path))
            return False

    else:
        logger.debug('performing a rmtree')
        rmtree(path)
        return True
