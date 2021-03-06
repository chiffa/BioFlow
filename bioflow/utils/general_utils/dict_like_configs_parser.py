"""
Tools for a more sane conversion between dicts and configfiles, based on the Python's default
ConfigParser
"""
import configparser


def ini_configs2dict(path):
    """
    Parses a config file on given path and converts it to a dict, in case of failure raises an IOError,

    :param path:
    :return:
    """
    stock_parser = configparser.SafeConfigParser()
    stock_parser_output = stock_parser.read(path)
    if not stock_parser_output:
        raise IOError('Cannot load configs file from %s' % path)
    parse_dict = {}
    for section in stock_parser.sections():
        parse_dict[section] = {}
        for option in stock_parser.options(section):
            parse_dict[section][option] = stock_parser.get(section, option)
    return parse_dict


def dict2init_configs(path, confdict):
    """
    Converts a dict of configs to a .ini files written to the location specified in the path argument
    :param path:
    :return:
    """
    base_writer = configparser.SafeConfigParser()
    for section, section_contents in confdict.items():
        base_writer.add_section(section)
        for parameter, parameter_value in section_contents.items():
            base_writer.set(section, parameter, parameter_value)

    with open(path, 'wt') as configfile:
        base_writer.write(configfile)
