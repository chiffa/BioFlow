__author__ = 'ank'

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from PolyPharma.configs import IDFilter, edge_type_filters

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


def Look_up_Annot_Node(p_load, p_type = ''):
    """
    return format: node type, node's displayName, node's db_ID, node's legacy ID

    .. code-block: python
    >>> print Look_up_Annot_Node('ENSG00000131981', 'UNIPROT_Ensembl')
    >>> # TODO: add the result here
    """

    def run_through(node_generator):
        """
        Gets the nodes that are referenced by the annot_nodes in the generator
        """
        retset = []
        for node in node_generator:
            gen_2 = node.inV("is_annotated")
            if not node_generator:
                raise Warning(str(node) + "is floating alone in the wild. He feels lonely.")
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
    pload = p_load.upper()
    if p_type == '':
        node_generator = DatabaseGraph.AnnotNode.index.lookup(payload = pload)
        if not node_generator:
            return []
        return run_through(node_generator)

    if p_type in Anot_Node_ptypes:
        node_generator =  double_index_search(pload, p_type)
        if not node_generator:
            return []
        return run_through(node_generator)

    raise Exception(p_type + "is unsupported. Please refer to Anot_Node_ptypes in neo4j_typeDec for supported types")


def Erase_custom_fields():
    """
        Resets the .costum field of all the Nodes on which we have iterated here. Usefull to perform
        after node set or node connectivity were modfied.

        Unlike the method in the Matrix_retrieval cluster, this method is very time-consuming, since it iterates
        on all the elements of all the classes susceptible to have the costum field.
    """

    # TODO: reconfigure to erase connexity infos later on.
    Node_gen = DatabaseGraph.Node.get(costum = 'Main_Connex')

    for Node in Node_gen:
        Node.custom = ''
        Node.save()


def reaction_participant_getter(Reaction, main_connex_only):
    edge_type_filter = edge_type_filters["Reaction"]
    LocalList = []
    count = 0
    for edge_type in edge_type_filter:
        if Reaction.bothV(edge_type) == None:
            continue
        for elt in Reaction.bothV(edge_type):
            Connex = True

            if main_connex_only:
                Connex = False
                if  elt.custom != None and "Main_Connex" in elt.custom:
                    Connex = True

            ID = str(elt).split('/')[-1][:-1]
            if ID not in IDFilter and Connex:
                LocalList.append(ID)
                count += 1
    return LocalList, count


def expand_from_seed(Seed_Node_ID, edge_filter):
    Seed_Node = DatabaseGraph.vertices.get(Seed_Node_ID)
    LocalList = []
    count = 0
    for edge_type in edge_filter:
        if Seed_Node.bothV(edge_type) != None:
            for elt in Seed_Node.bothV(edge_type):
                ID = str(elt).split('/')[-1][:-1]
                if ID not in IDFilter:
                    LocalList.append(ID)
                    count += 1
    return LocalList, count



if __name__ == "__main__":
    # print count_items(DatabaseGraph.UNIPORT)
    # lookup_by_ID(DatabaseGraph.UNIPORT, "CK2N2_HUMAN")
    # Erase_custom_fields()

    print Look_up_Annot_Node('ENSG00000131981', 'UNIPROT_Ensembl')