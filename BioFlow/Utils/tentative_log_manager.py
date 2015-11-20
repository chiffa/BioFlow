"""
File managing most of the logging behavior of the application.
"""
import os
import logging
import sys
from BioFlow.configs2 import log_location
from BioFlow.Utils.GeneralUtils.SanerFilesystem import mkdir_recursive

mkdir_recursive(log_location)  # create location where the code will be stored

logger = logging.getLogger('main_logger')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

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
