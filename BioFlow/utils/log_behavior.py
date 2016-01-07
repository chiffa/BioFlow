"""
File managing most of the logging behavior of the application.
"""
import os
from os import path
import logging
import sys
from shutil import rmtree


log_location = path.join(path.abspath(
    path.join(
        path.join(path.dirname(__file__), os.pardir),
        os.pardir)), 'logs')

on_unittest = os.environ.get('UNITTESTING') == 'True'  # if we are unittesting
on_remote_unittest = os.environ.get('REMOTE_UNITTEST') == 'True'  # if we are testing on CI tools
on_remote = os.environ.get('REMOTE') == 'True'

on_dev = False
# on_dev = True

# #################################
# redirecting all to a log file
# if on_remote = True:
#   f = open('../logs/Commander_logs.log','w')
#   sys.stdout = f
# ################################


def mkdir_recursive(my_path):  # pragma: no cover
    """
    Copy of mkdir recursive from saner configs, used here to remove circular dependencies
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param my_path:
    :return:
    """
    my_path = os.path.abspath(my_path)
    directory_name = os.path.dirname(my_path)
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(my_path):
        if '.' not in my_path.split(
                '/')[-1][-5:]:  # should be able to suppress specific file creation
            os.mkdir(my_path)


def wipe_dir(_path):  # pragma: no cover
    """
    wipes the indicated directory
    :param _path:
    :return: True on success
    """
    _path = os.path.abspath(_path)
    if not os.path.exists(_path):
        return True  # Nothing to do: destruction already done
    if os.path.isdir(_path):
        directory_name = _path
    else:
        directory_name = os.path.dirname(_path)
    if not os.path.isdir(directory_name):
        return False
    for sub_path in os.listdir(directory_name):
        if os.path.isdir(sub_path):
            return False
    rmtree(directory_name)
    return True


def add_handler(_logger, level, file_name):
    """
    Adds a file-writing handler for the log.

    :param _logger:
    :param level: logging.DEBUG or other level
    :param file_name: short file name, that will be stored within the application logs location
    :return:
    """
    _fh = logging.FileHandler(os.path.join(log_location, file_name), mode='a')
    _fh.setLevel(level)
    _fh.setFormatter(formatter)
    _logger.addHandler(_fh)


# define a formatter
formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

if on_dev:
    wipe_dir(log_location)

# create location where the code will be stored
mkdir_recursive(log_location)


def get_logger(logger_name):
    _logger = logging.getLogger(logger_name)
    _logger.setLevel(logging.DEBUG)

    add_handler(_logger, logging.DEBUG, 'debug.log')
    add_handler(_logger, logging.INFO, 'info.log')
    add_handler(_logger, logging.WARNING, 'warning.log')
    add_handler(_logger, logging.ERROR, 'error.log')
    add_handler(_logger, logging.CRITICAL, 'critical.log')

    if not on_remote_unittest:  # pragma: no cover
        ch = logging.StreamHandler(sys.stdout)
        if on_unittest or on_dev:
            ch.setLevel(logging.DEBUG)
        else:
            ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        _logger.addHandler(ch)

    return _logger

logger = get_logger('this_logger_needs_to_be_renamed')
