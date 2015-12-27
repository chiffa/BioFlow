"""
Module wrapping the class containing the repository retrieving engine
"""
from hashlib import md5, sha1
from BioFlow.utils.log_behavior import get_logger
from BioFlow.utils.general_utils.internet_io import url_to_local, check_hash

log = get_logger()


class DataRepository(object):
    """
    Single data repository object wrapper
    """

    def __init__(self, url, local_location,
                 checksum=None, checksum_type='md5', reference_doi=None):
        self.url = url
        self.local_location = local_location
        self.checksum = checksum
        if checksum_type in ['md5', 'sha1']:
            if checksum_type == 'md5':
                self.checksum_type = md5
            elif checksum_type == 'sha1':
                self.checksum_type = sha1
        else:
            log.warning('Cannot check file hash for %s from %s: checksum is not an md5 or sha1',
                        local_location, url)
            self.checksum = None

        # todo: generic organism names, organism selectors

    def pull(self):
        """
        Copies the file from online repository to
        """
        url_to_local(self.url, self.local_location)
        if self.checksum:
            checksum_matched = check_hash(self.local_location, self.checksum, self.checksum_type)
            if checksum_matched:
                log.info('successfully retrieved and checked %s from %s',
                         self.local_location,
                         self.url)
                return True
            else:
                log.info('Checksum mismatch for %s from %s',
                         self.local_location,
                         self.url)
                return False
        else:
            log.info('retrieved %s from %s; no checksum available',
                     self.local_location,
                     self.url)


class RepoEngine(object):
    """
    Parses a set of configuration rules and puts them all into repo wrapper objects
    """

    def __init__(self):
        pass