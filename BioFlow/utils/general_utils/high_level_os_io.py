"""
Saner io and filesystem manipulation compared to python defaults
"""
import os
from shutil import rmtree
from BioFlow.utils.log_behavior import logger


def mkdir_recursive(path):  # pragma: no cover
    """
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param path:
    :return:
    """
    logger.debug(
        'trying to create recursively path containing: %s', path)
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    logger.debug('subpath: %s', directory_name)
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(path):
        logger.debug('path %s does not exist yet, looks like file: %s' %
                     (path, '.' in path.split('/')[-1][-5:]))
        if '.' not in path.split('/')[-1][-5:]:
            # should be able to suppress specific file creation
            os.mkdir(path)
            logger.debug('; created')
        else:
            logger.debug('; creation skipped')


def wipe_dir(path):  # pragma: no cover
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

    logger.debug('going to wipe %s for path %s', directory_name, path)
    if not os.path.isdir(directory_name):
        logger.exception(
            'failed to delete %s: for path %s, not a dir', directory_name, path)
        return False

    for sub_path in os.listdir(directory_name):
        if os.path.isdir(sub_path):
            logger.exception(
                'failed to delete %s: for path %s, anti rm -rf flag', directory_name, path)
            return False

    else:
        logger.debug('performing a rmtree')
        rmtree(directory_name)
        return True
