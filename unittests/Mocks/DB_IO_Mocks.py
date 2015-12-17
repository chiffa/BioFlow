from mock import Mock as MagicMock
from BioFlow.utils.log_behavior import get_logger

log = get_logger('DB_IO_Mocks @ unittests')
log.debug('Mocking DB_IO module')


def look_up_annotation_set(supplied_list):
    return supplied_list, [(elt, '') for elt in supplied_list], ['' for _ in supplied_list]

