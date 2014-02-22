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
    i = 0
    for elt in Domain.get_all():
        i += 1
    return i

def Look_up_by_ID_for_a_set(Domain, ID_set):
    for ID in ID_set:
        print "scanning for:", ID
        lookup_by_ID(Domain, ID)
        print "=================="

# TODO: reset all the values to convert to the uppercase convetion.
def Look_up_Annot_Node(p_load, p_type=''):
    """
    return format: node type, node's displayName, node's db_ID, node's legacy ID
    """

    def run_through(node_generator):
        """
        Gets the nodes that are referenced by the annot_nodes in the generator
        """
        retset = []
        for node in node_generator:
            gen_2 = node.inV("is_annotated")
            if not node_generator:
                raise Warning(str(node)+"is floating alone in the wild. He feels lonely.")
            for rel_node in gen_2:
                node_db_ID = str(rel_node).split('/')[-1][:-1]
                node_ID = rel_node.ID
                node_type = rel_node.element_type
                node_display = rel_node.displayName
                retset.append((node_type, node_display, node_db_ID, node_ID))

        return retset

    def double_index_search(pload, ptype):
        node_generator = DatabaseGraph.AnnotNode.index.lookup(payload = pload)
        retset = []
        if node_generator:
            for node in node_generator:
                if ptype in node.ptype:
                    retset.append(node)
        return retset

    from PolyPharma.neo4j_Declarations.neo4j_typeDec import Anot_Node_ptypes
    p_load = p_load.upper()
    if p_type == '':
        node_generator = DatabaseGraph.AnnotNode.index.lookup(payload = p_load)
        if not node_generator:
            return []
        return run_through(node_generator)

    if p_type in Anot_Node_ptypes:
        node_generator =  double_index_search(p_load, p_type)
        if not node_generator:
            return []
        return run_through(node_generator)

    raise Exception(p_type + "is unsupported. Please refer to Anot_Node_ptypes in neo4j_typeDec for supported types")



if __name__ == "__main__":
    # print count_items(DatabaseGraph.UNIPORT)
    # lookup_by_ID(DatabaseGraph.UNIPORT,"CK2N2_HUMAN")

    print Look_up_Annot_Node('ENSG00000131981', 'UNIPROT_Ensembl')