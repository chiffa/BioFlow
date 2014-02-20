__author__ = 'ank'

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph

def lookup_by_ID(Domain, req):
    retset = Domain.index.lookup( ID = req )
    if not retset:
        print "nothing found"
    if retset:
        for item in retset:
            print item

def count_items(Domain):
    i=0
    for elt in Domain.get_all():
        i+=1
    return i

def Look_up_by_ID_for_a_set(Domain, ID_set):
    for ID in ID_set:
        print "scanning for:", ID
        lookup_by_ID(Domain, ID)
        print "=================="


Tamara_analist = ["CK2N2_HUMAN", "NPTX1_HUMAN", "PGS2_HUMAN", "CO8A1_HUMAN", "MALL_HUMAN", "EMAL1_HUMAN"]

if __name__ == "__main__":
    # print count_items(DatabaseGraph.UNIPORT)
    # lookup_by_ID(DatabaseGraph.UNIPORT,"CK2N2_HUMAN")
    Look_up_by_ID_for_a_set(DatabaseGraph.UNIPORT, Tamara_analist)
