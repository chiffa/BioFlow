"""

:author: Andrei Kucharavy
"""

__author__ = 'ank'

from PolyPharma.neo4j_Declarations.Graph_Declarator import DatabaseGraph
from PolyPharma.configs2 import IDFilter, Leg_ID_Filter, edge_type_filters, Dumps, Outputs, Hits_source, prename1, prename2
from PolyPharma.configs2 import Background_source, bgList
import pickle
from pprint import PrettyPrinter
from collections import defaultdict
from csv import reader, writer
import os


pp = PrettyPrinter(indent = 4)


def lookup_by_ID(Domain, req):
    """
    Looks up a node by legacy ID and prints it's access url. The lookup is a case-sensitive strict match.

    :param Domain: Node type of the form of DatabaseGraph.Object
    :param req: requested legacy ID
    """
    accumulator = []
    retset = Domain.index.lookup( ID = req )
    if not retset:
        print "nothing found"
    if retset:
        for item in retset:
            print item
            accumulator.append(item)


def count_items(Domain):
    """
    Stupid counter that gets the number of items in a given domain

    :warning: this method is highly underoptimal

    :param Domain: Domain of the DatabaseGraph.Object whose objects we are going to count
    :return: number of objects in the domain
    """
    i = 0
    for elt in Domain.get_all():
        i += 1
    return i


def Look_up_by_ID_for_a_set(Domain, ID_set):
    """
    Looks up nodes by legacy IDs from an ID_set list and prints their access url. The lookup is a case-sensitive strict match.

    :param Domain: Node type of the form of DatabaseGraph.Object
    :param Domain: Node type of the form of DatabaseGraph.Object
    :param ID_set: list of requested legacy ID
    """
    for ID in ID_set:
        print "scanning for:", ID
        lookup_by_ID(Domain, ID)
        print "=================="


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
    :raise Warning: if an annotation node is not bound to a real node. This might happen if some object nodes were manually
                    deleted, but not their annotation nodes. Tu curb this a full database reload is required
    """
    if not node_generator:
        return []

    retset = []
    for node in node_generator:
        gen_2 = node.inV("is_annotated")
        if not gen_2:
            raise Warning(str(node) + "is floating alone in the wild. He feels lonely.")
        for rel_node in gen_2:
            node_db_ID = str(rel_node).split('/')[-1][:-1]
            node_ID = rel_node.ID
            node_type = rel_node.element_type
            node_display = rel_node.displayName
            retset.append((node_type, node_display, node_db_ID, node_ID))
    return retset


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



def look_up_Annot_Node(p_load, p_type = ''):
    """
    Looks up nodes accessible via the annotation nodes with a given annotation and given annotation type.
    The lookup strict match, but case-insensitiYOR031Wve.

    .. code-block: python
    >>> print look_up_Annot_Node('ENSG00000131981', 'UNIPROT_Ensembl')
    >>> # TODO: add the results here when implementing the doctests


    :param p_load: payload
    :param p_type: payload type
    :return: node type, node's displayName, node's db_ID, node's legacy ID
    :rtype: 4- tuple
    :raise Exception: "p_type unsupported", in case a p_type is not on the supported list specified in the neo4j_Declarations.neo4j_typeDec
    """

    def double_index_search(pload, ptype):
        """
        Supporting fucntion. Performs a search in the annotation nodes over both a payload and a ptype

        :param pload: payload
        :param ptype: payload type
        :return: list of found annotation nodes satisfying both payload content and payload type conditions
        """
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
        return run_through(node_generator)

    if p_type in Anot_Node_ptypes:
        node_generator =  double_index_search(pload, p_type)
        return run_through(node_generator)

    raise Exception(p_type + "is unsupported. Please refer to Anot_Node_ptypes in neo4j_typeDec for supported types")


def look_up_Annot_set(p_load_list, p_type=''):
    """
    Looks up an set of annotations in the database and finds the Ids of nodes containing SWISSPROT proteins

    :param p_load_list:
    :param p_type:
    :return:
    """
    retdict = dict( (p_load, look_up_Annot_Node(p_load, p_type)) for p_load in p_load_list)
    retlist = [value[0][2] for value in retdict.itervalues() if value!=[]]
    warnlist = [key for key, value in retdict.iteritems() if value == []]
    for warnId in warnlist:
        print Warning('following ID has no correspondance in the database: ' + warnId)
    return warnlist, retdict, retlist


def Erase_custom_fields():
    """
        Resets the .costum field of all the Nodes on which we have iterated here. Usefull to perform
        after node set or node connectivity were modfied.
        Unlike the method in the Matrix_retrieval cluster, this method is very time-consuming, since it iterates
        on all the elements of all the classes susceptible to have the costum field.
    """
    Node_gen = DatabaseGraph.vertices.index.lookup(costum = 'Main_Connex')
    if Node_gen:
        for Node in Node_gen:
            Node.custom = ''
            Node.main_connex = False
            Node.save()

    Node_gen = DatabaseGraph.vertices.index.lookup(main_connex = True)
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
        if Reaction.bothV(edge_type) == None:
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
        if Seed_Node.bothV(edge_type) != None:
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
            generator = bulbs_type.index.lookup(displayName = forbidden_Legacy_ID)
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
        generator = DatabaseGraph.UNIPORT.index.lookup(ID = node_leg_Id)
        if not generator:
            continue
        retlist = []
        for node in generator:
            retlist.append(node.displayName)
        if len(retlist)!=1:
            raise Exception('Something went wrong with the UP retrieval for the UP %s, too many display names: %s' %
                            (node_leg_Id,retlist))
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
    retlist=[]
    with open(Hits_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row
    retlist = [ret for ret in retlist]
    source = look_up_Annot_set(retlist)
    PrettyPrinter(indent=4, stream=open(prename1, 'w')).pprint(source[1])
    writer(open(prename2, 'w'), delimiter='\n').writerow(source[2])


def unwrap_background():
    retlist=[]

    with open(Background_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row

    with open(Hits_source) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            retlist = retlist + row

    retlist = list(set(ret for ret in retlist))
    source = look_up_Annot_set(retlist)
    writer(open(bgList, 'w'), delimiter='\n').writerow(source[2])


# TODO: should be refactored into the configs
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    Forbidden_verification_dict = {   'Small Molecule': DatabaseGraph.SmallMolecule,
                                  'Small Molecule Collection': DatabaseGraph.SmallMolecule_Collection,
                                  'Physical Entity': DatabaseGraph.PhysicalEntity,
                                  'Physical Entity Collection': DatabaseGraph.PhysicalEntity_Collection,
                                }
else:
    Forbidden_verification_dict = { }

if __name__ == "__main__":

    transcription = ['1005777', '1001874', '842142', '836106', '1014143', '1058552', '821021', '826066', '886586',
                     '865200', '835714', '766912', '900540', '811360', '816278', '1079377', '740496', '1002986',
                     '778615', '1000616', '950906', '996907', '1033828', '971819', '809996', '781060', '874882', '758558',
                     '740148', '758263', '926744', '907135', '990196', '743186', '1011403', '979456', '1073345', '913095',
                     '805935', '777825', '1028447', '951846', '1075746', '882269', '888669', '1029555', '941071', '751245',
                     '986669', '927801', '783375', '1066551', '797230', '859473', '914620', '1083041', '837861', '901826',
                     '906807', '913721', '991031', '831309', '810280', '856443', '986631', '800198', '832809', '774804',
                     '1048252', '935385', '92{0480', '952510', '988829', '744318', '1001042', '1069894', '848924', '949209',
                     '1035740', '770322', '749382', '1003600', '889200', '1071170', '841076', '822788', '810700', '801219',
                     '801191', '964332', '935139', '1071008', '1002847', '1022239', '875186', '848997', '966988', '1008725',
                     '899845', '763208', '823114', '1078022', '969713', '892913', '1104135', '794981', '1005373', '767473',
                     '820454', '1066442', '753200', '1000686', '800665', '991315', '817329', '885986', '1019819', '871354',
                     '1073014', '886044', '1004612', '871049', '1034170', '860442', '848198', '878608', '844231', '794289',
                     '1016545', '1023218', '1103896', '785577', '995803', '912780', '897083', '886840', '750127', '945831',
                     '916664', '864667', '746441', '1039958', '1101959', '956599', '875904', '1073318', '909015', '974199',
                     '1037879', '964653', '784732', '949392', '1091093', '1081714', '837765', '973051', '772194', '765760',
                     '973374', '893485', '838733', '942660', '798235', '996784', '914154', '1055665', '935238', '780699',
                     '1016715', '991271', '955633', '895017', '794788', '1091063', '1046353', '986976', '826048', '859213',
                     '839038', '783622', '760894', '853760', '1106249', '819917', '1062291', '815494', '862693', '814921',
                     '938485', '783138', '885381', '903943', '772264', '980817', '823252', '962633', '758463', '755526',
                     '757978', '765208', '863431', '1042302', '905303', '811490', '939476', '846291', '993567', '819257',
                     '886290', '1013510', '827509', '991517', '1026799', '766799', '1029295', '806975', '911764', '1070622',
                     '891117', '782645', '994623', '1013051', '830556', '844362', '1040267', '933951', '876062', '1060754',
                     '826120', '973468', '949926', '789918', '747228', '764653', '1056009', '893389', '828460', '780732',
                     '1061892', '748295', '890461', '800707', '896612', '1061111', '1038976', '1026861', '892596', '994512',
                     '908624', '819172', '949891', '1015019', '802349', '752772', '1041116', '920885', '851881', '758797',
                     '898171', '844767', '775297', '936117', '1089489', '970498', '854601', '932389', '911677', '773272',
                     '964086', '837996', '968860', '954895', '847543']


    translation = ['860713', '1054645', '1050811', '1010442', '746034', '960074', '841626', '989237', '992333', '885431',
                   '1049088', '1005233', '742678', '1025641', '1081901', '940598', '1004992', '738650', '1018277',
                   '997379', '745984', '1051693', '758243', '1085622', '933193', '1070658', '932056', '774142',
                   '1001023', '1030061', '743020', '905602', '998094', '932332', '996254', '823270', '856693',
                   '1017722', '976875', '960939', '890389', '1099836', '829293', '957345', '1020969', '905383',
                   '1009002', '771146', '1081329', '985624', '904324', '1019095', '937024', '1051946', '854456',
                   '1076560', '741935', '1102657', '1021881', '1041371', '757039', '997830', '942299', '773580',
                   '902901', '1069785', '1023488', '1077664', '864199', '815792', '947643', '983830', '1049137',
                   '821120', '1094403', '958096', '877347', '870065', '779276', '998107', '966557', '876290', '782349',
                   '886637', '973293', '943484', '840896']



    GBO_1 = ['583954', '565151', '625184', '532448', '553020', '547608', '576300', '533299', '540532', '591419']
    #     [
    # 'YOR031W',
    # 'YOR001W',
    # 'YOL107W',
    # 'YOL124C',
    # 'YOL040C',
    # 'YOR184W',
    # 'YOR374W',
    # 'YOR125C',
    # 'YOL087C',
    # 'YOR334W',
    # ]
    
    GBO_2 = ['562293', '544722', '534354', '612635', '532463', '561658', '630018', '586185', '611762', '599295']
    #     [
    # 'YOL084W',
    # 'YOR243C',
    # 'YOR281C',
    # 'YOR127W',
    # 'YOR250C',
    # 'YOR278W',
    # 'YOR354C',
    # 'YOR319W',
    # 'YOL114C',
    # 'YOR271C',
    # ]
    
    GBO_3 = [
    'YOL040C',
    'YOL124C',
    'YOL155C',
    'YOR031W',
    'YOR080W',
    'YOR253W',
    'YOR316C-A',
    'YOR332W',
    'YOR338W',
    'YOR339C',
    ]
    
    GBO_4 = ['594353', '565151', '618791', '537788', '546413', '576300', '533299', '540532', '532448', '557819']
    #     [
    # 'YOR339C',
    # 'YOR316C-A',
    # 'YOR184W',
    # 'YOR167C',
    # 'YOR138C',
    # 'YOR031W',
    # 'YOR001W',
    # 'YOL124C',
    # 'YOL040C',
    # 'YOL015W',
    # ]


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
    pass