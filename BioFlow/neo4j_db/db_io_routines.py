"""
Module responsible for the declaration of all the database access routines that are used by other
modules. In case a different back-end than neo4j is used, all methods in this cluster have to be
re-implemented
"""
import os
import pickle
from collections import defaultdict
from csv import reader, writer
from pprint import PrettyPrinter

from BioFlow.neo4j_db.GraphDeclarator import DatabaseGraph
from BioFlow.main_configs import IDFilter, Leg_ID_Filter, edge_type_filters, Dumps, Hits_source, \
    prename1, prename2, Background_source, bgList, annotation_nodes_ptypes

pp = PrettyPrinter(indent=4)


def lookup_by_id(domain, req):
    """
    Looks up a node by legacy ID and prints it's access url. The lookup is a case-sensitive strict
    match.

    :param domain: Node type of the form of DatabaseGraph.Object
    :param req: requested legacy ID
    :return: list of the items found
    """
    accumulator = []
    retset = domain.index.lookup(ID=req)
    if not retset:
        print "nothing found"
    if retset:
        for item in retset:
            print item
            accumulator.append(item)

    return accumulator


def count_items(domain):
    """
    Stupid counter that gets the number of items in a given domain

    :param domain: Domain of the DatabaseGraph.Object whose objects we are going to count
    :return: number of objects in the domain
    """
    return sum(1 for _ in domain.get_all())


def set_look_up_by_id(domain, id_set):
    """
    Looks up nodes by legacy IDs from an ID_set list and prints their access url. The lookup is a
    case-sensitive strict match.

    :param domain: Node type of the form of DatabaseGraph.Object
    :param Domain: Node type of the form of DatabaseGraph.Object
    :param id_set: list of requested legacy ID
    """
    for ID in id_set:
        print "scanning for:", ID
        lookup_by_id(domain, ID)
        print "=================="

# Looks like the two above nodes are


def get_attached_annotations(node_id):
    retset = []
    node = DatabaseGraph.vertices.get(node_id)
    gen_2 = node.outV("is_annotated")
    if not gen_2:
        print "node %s has not annotations attached to it" % node_id
    for rel_node in gen_2:
        node_db_ID = str(rel_node).split('/')[-1][:-1]
        retset.append(node_db_ID)
    return retset


def run_through(node_generator):
    """
    Supporting function. Gets the nodes that are referenced by the annot_nodes in the generator

    :param node_generator: iterator over annot_nodes
    :return: the list of object nodes accessible from this set of annot_nodes
    """

    if not node_generator:  # Node generator is empty: return empty list
        return []

    set_of_interest = []
    for node in node_generator:
        gen_2 = node.inV("is_annotated")
        # backlink to the node which this node annotates
        if not gen_2:
            # if not backlink was found, raise a warning about lonely node
            print Warning("%s is floating alone in the wild. He feels lonely." % str(node))
        for rel_node in gen_2:
            node_db_id = str(rel_node).split('/')[-1][:-1]
            node_id = rel_node.ID
            node_type = rel_node.element_type
            node_display = rel_node.displayName
            set_of_interest.append(
                (node_type, node_display, node_db_id, node_id))
    return set_of_interest


def unwrap_DB_ID(node_generator):
    """
    Supporting function. Gets the nodes that are referenced by the annot_nodes in the generator

    :param node_generator: iterator over annot_nodes
    :return: the DB_IDs list of object nodes accessible from this set of annot_nodes
    :raise Warning: if an annotation node is not bound to a real node. This might happen if some object nodes were manually
                    deleted, but not their annotation nodes. Tu curb this a full database reload is required
    """
    if not node_generator:
        return []

    retset = []
    for node in node_generator:
        node_db_ID = str(node).split('/')[-1][:-1]
        retset.append(node_db_ID)
    return retset


def look_up_annotation_node(p_load, p_type=''):
    """
    Looks up nodes accessible via the annotation nodes with a given annotation and given
    annotation  type.
    The lookup strict match, but case-insensitive.

    :param p_load: payload
    :param p_type: payload type. servers to restrict search to a specific ID subset
    :return: (node type, node's displayName, node's db_ID, node's legacy ID)
    :raise Exception: "p_type unsupported", in case a p_type is not on the supported list
    specified  in the neo4j_db.neo4j_typeDec
    """

    def double_index_search(pay_load, payload_type):
        """
        Supporting function. Performs a search in the annotation nodes over both a payload and a
        payload type

        :param pay_load: payload
        :param payload_type: payload type
        :return: list of found annotation nodes satisfying both payload content and payload type
         conditions
        """
        inner_node_generator = DatabaseGraph.AnnotNode.index.lookup(
            payload=pay_load)
        results = []
        if inner_node_generator:  # if there is anything in the inner node generator
            for node in inner_node_generator:
                if payload_type in node.ptype:  # late filtering by the node ptype
                    results.append(node)
        return results

    payload = p_load.upper()  # make payload case-agnostic

    if p_type in annotation_nodes_ptypes:
        node_generator = double_index_search(payload, p_type)
        return run_through(node_generator)

    if p_type != '':
        print Exception("%s is an unsupported payload type. Please refer to " +
                        "annotation_nodes_ptypes in cofigs/internal_configs.py. ignoring p_type " +
                        "from now on")

    node_generator = DatabaseGraph.AnnotNode.index.lookup(payload=payload)
    return run_through(node_generator)


def look_up_annotation_set(p_load_list, p_type=''):
    """
    Looks up an set of annotations in the database and finds the Ids of nodes containing SWISSPROT
    proteins linked to by annotations

    :param p_load_list: list of payloads
    :param p_type: expected type of payloads
    :return:
    """
    def db_id_mapping_helper(mapped_db_id_list):
        if mapped_db_id_list:
            return mapped_db_id_list[0][2]
        else:
            return ''

    tuple_list = [(p_load, look_up_annotation_node(p_load, p_type))
                  for p_load in p_load_list]
    db_id_list = [db_id_mapping_helper(value) for key, value in tuple_list]
    warnings_list = [key for key, value in tuple_list if value == []]
    for warnId in warnings_list:
        print Warning('look_up_annotation_set@DB_IO: following ID has no corresponding entry in ' +
                      'the database: %s' % warnId)

    # TODO:
    return warnings_list, tuple_list, db_id_list


def erase_custom_fields():
    """
        Resets the .costum field of all the Nodes on which we have iterated here. Usefull to perform
        after node set or node connectivity were modfied.
        Unlike the method in the Matrix_retrieval cluster, this method is very time-consuming, since it iterates
        on all the elements of all the classes susceptible to have the costum field.
    """
    Node_gen = DatabaseGraph.vertices.index.lookup(costum='Main_Connex')
    if Node_gen:
        for Node in Node_gen:
            Node.custom = ''
            Node.main_connex = False
            Node.save()

    Node_gen = DatabaseGraph.vertices.index.lookup(main_connex=True)
    if Node_gen:
        for Node in Node_gen:
            Node.custom = ''
            Node.main_connex = False
            Node.save()


def reaction_participant_getter(Reaction, main_connex_only):
    """
    Recovers all the participants of the reaction

    :param Reaction: Reaction node for which we are willing to get the participants
    :param main_connex_only: If set to true, will only pull elements from the reaction that are in the main connex_set
    :type main_connex_only: bool
    :return: List of found nodes, number of found nodes
    :rtype: 2- tuple
    """
    edge_type_filter = edge_type_filters["Reaction"]
    LocalList = []
    count = 0
    for edge_type in edge_type_filter:
        if Reaction.bothV(edge_type) is None:
            continue
        for elt in Reaction.bothV(edge_type):
            Connex = True

            if main_connex_only:
                Connex = False
                if elt.main_connex:
                    Connex = True

            ID = str(elt).split('/')[-1][:-1]
            if ID not in IDFilter and Connex:
                LocalList.append(ID)
                count += 1
    return LocalList, count


def expand_from_seed(Seed_Node_ID, edge_filter, main_connex_only):
    """
    Recovers all the nodes accessible in one jump from a seed_node with a given database ID by jumping only via the relations
        of type specified in the edge_filter

    :param Seed_Node_ID: the database ID of the initial node from which we are observing accessibility
    :param edge_filter: the list of relation types for which the jumps are authorised
    :return: List of found nodes database ID, number of found nodes
    :rtype: 2- tuple
    """
    Seed_Node = DatabaseGraph.vertices.get(Seed_Node_ID)
    LocalList = []
    count = 0
    for edge_type in edge_filter:
        if Seed_Node.bothV(edge_type) is not None:
            for elt in Seed_Node.bothV(edge_type):
                Connex = True

                if main_connex_only:
                    Connex = False
                    if elt.main_connex:
                        Connex = True

                ID = str(elt).split('/')[-1][:-1]
                if ID not in IDFilter and Connex:
                    LocalList.append(ID)
                    count += 1
    return LocalList, count


def recompute_forbidden_IDs(Node_Type_Dict):
    """
    Recomputes the list of nodes that contain overloaded terms that would bring too close together the reactions that are normally not,
    just because of participation of ultra-aboundant elements, such as H2O, H+ or ATP

    :param Node_Type_List: Dictionary mapping the names of entities to their corresponding bulbs classes
    :type Node_Type_Dict: dict
    """
    retlist = set()
    for bulbs_type in Node_Type_Dict.itervalues():
        for forbidden_Legacy_ID in Leg_ID_Filter:
            generator = bulbs_type.index.lookup(
                displayName=forbidden_Legacy_ID)
            UNW = unwrap_DB_ID(generator)
            retlist.update(UNW)
    print retlist
    print Dumps.Forbidden_IDs
    pickle.dump(retlist, file(Dumps.Forbidden_IDs, 'w'))


def recover_UP_chars(UP_Nodes, UP_are_IDs):
    retdict = {}

    if UP_are_IDs:
        for node_Id in UP_Nodes:
            node = DatabaseGraph.UNIPORT.get(node_Id)
            retdict[node] = [node.ID, node.displayName]
        return retdict

    for node_leg_Id in UP_Nodes:
        generator = DatabaseGraph.UNIPORT.index.lookup(ID=node_leg_Id)
        if not generator:
            continue
        retlist = []
        for node in generator:
            retlist.append(node.displayName)
        if len(retlist) != 1:
            raise Exception(
                'Something went wrong with the UP retrieval for the UP %s, too many display names: %s' %
                (node_leg_Id, retlist))
        else:
            retdict[node_leg_Id] = retlist
    return retdict


def recover_annotation(Node_Id_set, annotation_type):

    ret_dict = defaultdict(list)
    ret_list = []

    for node_id in Node_Id_set:
        node = DatabaseGraph.vertices.get(node_id)
        if node:
            annotation_generator = node.bothV('is_annotated')
            if annotation_generator:
                for annot in annotation_generator:
                    if annot.ptype == annotation_type:
                        if 'G0' in annot.payload:
                            ret_dict[node_id].append(annot.payload)
                            ret_list.append([str(annot.payload), str(node.ID)])

    return dict(ret_dict), ret_list


def unwrap_source():
    retlist = []
    with open(Hits_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row
    retlist = [ret for ret in retlist]
    source = look_up_annotation_set(retlist)
    PrettyPrinter(indent=4, stream=open(prename1, 'w')).pprint(source[1])
    writer(open(prename2, 'w'), delimiter='\n').writerow(source[2])


def unwrap_background():
    retlist = []

    with open(Background_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row

    with open(Hits_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row

    retlist = list(set(ret for ret in retlist))
    source = look_up_annotation_set(retlist)
    writer(open(bgList, 'w'), delimiter='\n').writerow(source[2])

# Yes, I know what goes below here is ugly and shouldn't be in the
# production part of the code

on_rtd = os.environ.get('READTHEDOCS') == 'True'
on_unittest = os.environ.get('UNITTESTING') == 'True'

if not on_rtd:
    Forbidden_verification_dict = {
        'Small Molecule': DatabaseGraph.SmallMolecule,
        'Small Molecule Collection': DatabaseGraph.SmallMolecule_Collection,
        'Physical Entity': DatabaseGraph.PhysicalEntity,
        'Physical Entity Collection': DatabaseGraph.PhysicalEntity_Collection, }
else:
    Forbidden_verification_dict = {}

if on_unittest:
    # Yes, this is dangerous as hell. Can't see a better way of doing it
    # though for now.
    import sys
    import unittests.Mocks.DB_IO_Mocks as self_mock
    sys.modules[__name__] = self_mock

# TODO: why are we manipulating the Forbidden_verification dictionary and Reactome import
# dictionary? Is it because they are used elsewhere and we need to cause a
# skip in code elsewhere?


if __name__ == "__main__":

    # print count_items(DatabaseGraph.UNIPORT)
    # lookup_by_ID(DatabaseGraph.UNIPORT, "SIR2_YEAST")
    # Erase_custom_fields()
    # recompute_forbidden_IDs(Forbidden_verification_dict)
    # print unwrap_source()
    # print look_up_Annot_Node('CTR86')
    # print Look_up_Annot_Node('ENSG00000131981', 'UNIPROT_Ensembl')
    # unwrap_source()
    # unwrap_background()
    # anset = [GBO_1, GBO_2, GBO_3, GBO_4]

    # for subset in anset:
    #     print look_up_Annot_set(subset)[-1]

    # print len(transcription)
    # print recover_annotation(transcription, 'UNIPROT_Ensembl')[1]
    # print recover_UP_chars(UP_Nodes=transcription, UP_are_IDs=None)
    # _, resdict, reslist = look_up_Annot_set(['MYPN'])
    # pp.pprint(resdict)
    # print reslist
    unwrap_source()
