"""
File managing most of the logging behavior of the application.
"""
import os
# from os import path
import logging
import logging.handlers
import sys
from shutil import rmtree
import traceback
import warnings
import threading

from bioflow.configs.bioflow_home import output_location, log_location

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


def _warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    """
    Supporting function to allow warnings with traceback

    :param message:
    :param category:
    :param filename:
    :param lineno:
    :param file:
    :param line:
    :return:
    """

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = _warn_with_traceback


def mkdir_recursive(my_path):  # pragma: no cover
    """
    Copy of mkdir recursive from saner configs, used here to remove circular dependencies
    Recursively creates a directory that would contain a file given a windows-like filename
    (xxx.xxx) or a directory name

    :param my_path: path which to recursively create
    :return: None
    """
    my_path = os.path.abspath(my_path)
    directory_name = os.path.dirname(my_path)
    if not os.path.exists(directory_name):
        mkdir_recursive(directory_name)
    if not os.path.exists(my_path):
        if '.' not in my_path.split(
                '/')[-1][-5:]:  # should be able to suppress specific file creation
            os.mkdir(my_path)


def wipe_dir(my_path):  # pragma: no cover
    """
    wipes the indicated directory

    :param my_path: path to wipe
    :return: True on success
    """
    my_path = os.path.abspath(my_path)
    if not os.path.exists(my_path):
        return True  # Nothing to do: destruction already done
    if os.path.isdir(my_path):
        directory_name = my_path
    else:
        directory_name = os.path.dirname(my_path)
    if not os.path.isdir(directory_name):
        return False
    for sub_path in os.listdir(directory_name):
        if os.path.isdir(sub_path):
            return False
    rmtree(directory_name, ignore_errors=True)
    return True


def add_to_file_handler(my_logger, level, file_name, rotating=False, log_location=log_location):
    """
    Adds a file-writing handler for the log.

    :param my_logger: the logger to which add a handler
    :param level: logging.DEBUG or other level
    :param file_name: short file name, that will be stored within the application logs location
    :param rotating: if true, rotating file handler will be added.
    :return:
    """
    handler_name = os.path.join(log_location, file_name)
    if rotating:
        _fh = logging.handlers.RotatingFileHandler(handler_name, maxBytes=1e7, backupCount=3)
    else:
        _fh = logging.FileHandler(handler_name, mode='a')
    _fh.setLevel(level)
    _fh.setFormatter(formatter)
    my_logger.addHandler(_fh)


# define a formatter
formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s')

if on_dev:
    wipe_dir(log_location)


# create location where the logs will be stored
mkdir_recursive(log_location)
mkdir_recursive(output_location)


def get_logger(logger_name):
    """
    Returns a properly configured logger object

    :param logger_name: name of the logger object
    """
    _logger = logging.getLogger(logger_name)
    _logger.setLevel(logging.DEBUG)

    add_to_file_handler(_logger, logging.DEBUG, 'debug.log', rotating=True)
    add_to_file_handler(_logger, logging.INFO, 'info.log')
    add_to_file_handler(_logger, logging.INFO, 'run.log', log_location=output_location)
    add_to_file_handler(_logger, logging.WARNING, 'warning.log')
    add_to_file_handler(_logger, logging.ERROR, 'error.log')
    add_to_file_handler(_logger, logging.CRITICAL, 'critical.log')
    # add_to_file_handler()

    def handle_exception(exc_type, exc_value, exc_traceback):

        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return

        _logger.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
        # raise RuntimeError("Terminating the Exception handling")

    sys.excepthook = handle_exception

    def handle_thread_exception(exc_type, exc_value, exc_traceback):

        _logger.critical("Uncaught thread exception", exc_info=(exc_type, exc_value, exc_traceback))
        raise RuntimeError("Terminating the Exception handling")

    threading.excepthook = handle_thread_exception

    logging.captureWarnings(True)
    default_warn_logger = logging.getLogger("py.warnings")
    add_to_file_handler(default_warn_logger, logging.WARNING, 'dependencies_warnings.log')

    if not on_remote_unittest:  # pragma: no cover
        ch = logging.StreamHandler(sys.stderr)
        if on_unittest or on_dev:
            ch.setLevel(logging.DEBUG)
        else:
            ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        _logger.addHandler(ch)

    return _logger


def clear_logs():
    """
    Wipes all logs
    """
    wipe_dir(log_location)
