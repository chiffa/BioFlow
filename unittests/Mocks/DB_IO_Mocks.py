from mock import Mock as MagicMock
from BioFlow.Utils.LogManager import logger

logger.debug('Mocking DB_IO module')


def look_up_annotation_set(supplied_list):
    return supplied_list, [(elt, '') for elt in supplied_list], ['' for _ in supplied_list]

