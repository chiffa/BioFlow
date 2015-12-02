"""
File managing most of the logging behavior of the application.
"""
import os
from os import path
import logging
import sys

log_location = path.join(path.abspath(
    path.join(
        path.join(path.dirname(__file__), os.pardir),
        os.pardir)), 'logs')


def mkdir_recursive(path):
    """
    Copy of mkdir recursive from saner configs, used here to remove circular dependencies
    Recursively creates a directory that would contain a file given win-like filename (xxx.xxx)
    or directory name
    :param path:
    :return:
    """
    path = os.path.abspath(path)
    directory_name = os.path.dirname(path)
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(path):
        if '.' not in path.split(
                '/')[-1][-5:]:  # should be able to suppress specific file creation
            os.mkdir(path)


mkdir_recursive(log_location)  # create location where the code will be stored

logger = logging.getLogger('main_logger')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(os.path.join(log_location, 'debug.log'), mode='a')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

fh = logging.FileHandler(os.path.join(log_location, 'info.log'), mode='a')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

fh = logging.FileHandler(os.path.join(log_location, 'warning.log'), mode='a')
fh.setLevel(logging.WARNING)
fh.setFormatter(formatter)
logger.addHandler(fh)

fh = logging.FileHandler(os.path.join(log_location, 'error.log'), mode='a')
fh.setLevel(logging.ERROR)
fh.setFormatter(formatter)
logger.addHandler(fh)

fh = logging.FileHandler(os.path.join(log_location, 'critical.log'), mode='a')
fh.setLevel(logging.CRITICAL)
fh.setFormatter(formatter)
logger.addHandler(fh)
