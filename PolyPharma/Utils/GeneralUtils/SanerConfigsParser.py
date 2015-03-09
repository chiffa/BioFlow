__author__ = 'ank'

import ConfigParser

def ini_configs2dict(path,):
    '''
    Parses a config file on given path, in case of failure raises an IOError,
    '''
    base_parser = ConfigParser.SafeConfigParser()
    parsed = base_parser.read(path)
    if parsed == []:
        raise IOError('Cannot load configs file from %s' % path)
    parse_dict = {}
    for section in base_parser.sections():
        parse_dict[section] = {}
        for option in base_parser.options(section):
            parse_dict[section][option] = base_parser.get(section, option)
    return parse_dict