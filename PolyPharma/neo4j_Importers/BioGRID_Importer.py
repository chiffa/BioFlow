__author__ = 'ank'

from csv import reader
from PolyPharma.configs2 import BioGRID
from PolyPharma.neo4j_analyzer.DB_IO_Routines import look_up_Annot_set
from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph


def parse_BioGRID(_BioGRID):
    """
    Parses the given file as a BioGrid file and returns as

    :param _BioGRID:
    :return:
    """
    ret_dict = {}
    base = []

    with open(_BioGRID,'rb') as sourcefile:
        rdr = reader(sourcefile, 'excel-tab')
        rdr.next()
        for row in rdr:
            ret_dict[tuple(row[7:9])] = [row[17]]
            if row[18] != '-':
                ret_dict[tuple(row[7:9])].append(row[18])
            base.append(row[7])
            base.append(row[8])

    return ret_dict, base


def convert_to_DB_Ids(base):
    warnlist, resdict, reslist = look_up_Annot_set(set(base))
    retdict = dict((key, value[0][2]) for key, value in resdict.iteritems() if key not in warnlist)
    print 'BioGrid ID cast converter length: %s' % len(retdict)
    return retdict


def cast_into_DB(re_caster, ret_dict):
    final_dicts = dict(((re_caster[key[0]],re_caster[key[1]]),value) for key, value in ret_dict.iteritems()
                       if key[0] in re_caster.keys() and key[1] in re_caster.keys())

    for (node_name1, node_name2), arglist in final_dicts.iteritems():
        node1 = DatabaseGraph.UNIPORT.get(node_name1)
        node2 = DatabaseGraph.UNIPORT.get(node_name2)
        if len(arglist) > 1:
            DatabaseGraph.is_weakly_interacting.create(node1,node2, throughput=arglist[0], confidence=float(arglist[1]))
        else:
            DatabaseGraph.is_weakly_interacting.create(node1,node2, throughput=arglist[0])


def import_BioGRID():
    ret_dict, base = parse_BioGRID(BioGRID)
    re_caster = convert_to_DB_Ids(base)
    cast_into_DB(re_caster, ret_dict)


if __name__ == "__main__":
    pass