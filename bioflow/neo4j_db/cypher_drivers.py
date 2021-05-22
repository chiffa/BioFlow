"""
Where the actual heavy-lifting of teh database connection is actually done.
Using a different database would require re-implementing the GraphDBPipe class to fit the
signatures.
"""
import os
from pprint import pprint
from collections import defaultdict
from itertools import combinations_with_replacement
from neo4j import GraphDatabase, DEFAULT_DATABASE
from neo4j.graph import Node, Relationship, Path
from typing import List, Tuple, NewType, Dict
from bioflow.utils.log_behavior import get_logger
from bioflow.configs.main_configs import neo4j_server_url, neo4j_db_name


log = get_logger(__name__)

# default connection parameters
uri = neo4j_server_url
user = 'neo4j'

# Type hinting support
db_id = NewType('db_id', int)
db_n_type = NewType('db_n_type', str)
# REFACTOR: [typing] add a db_a_type for annotation node and db_p_e_type for physical entities?
db_e_type = NewType('db_e_type', str)
e_orientation = NewType('e_orientation', str)


# This hard-coded, because those are essential to the logic of the application and traversal of
# the neo4j graph when parsing to laplace conversion
allowed_node_parse_types = ['physical_entity', 'annotation', 'xref']
allowed_edge_parse_types = ['physical_entity_molecular_interaction', 'identity', 'refines',
                            'annotates', 'annotation_relationship', 'xref']

required_node_params = {'parse_type', 'source', 'legacyID', 'displayName'}
required_edge_params = {'parse_type', 'source'}
# annotation tag does not get checked, because it is inserted by a separate function than the
# ones that does node/edge creation


def _check_node_params(param_dict: dict) -> bool:
    """
    Auxilary function

    Checks if the parameters provided to create a new node respect the constraints on the allowed
    `parse_tyoe`s and provide the required parameters

    :param param_dict:
    :return: True if all is good
    :raise Exception: exception in case either of the preconditions fail
    """
    if param_dict is None:
        raise Exception('param_dict supplied is empty')

    if required_node_params.issubset(param_dict.keys()):
        if param_dict['parse_type'] in allowed_node_parse_types:
            return True
        else:
            raise Exception('param_dict has a non-allowed parse_type: %s. Allowed: %s'
                            % (param_dict['parse_type'], allowed_node_parse_types))
    else:
        raise Exception('Required parameter(s) missing : %s'
                        % (required_node_params - param_dict.keys()))


def _check_edge_params(param_dict: dict) -> bool:
    """
    Auxilary function

    Checks if the parameters provided to create a new edge respect the constraints on the allowed
    `parse_tyoe`s and provide the required parameters

    :param param_dict:
    :return: True if all is good
    :raise Exception: exception in case either of the preconditions fail
    """
    if param_dict is None:
        raise Exception('param_dict supplied is empty')

    if required_edge_params.issubset(param_dict.keys()):
        if param_dict['parse_type'] in allowed_edge_parse_types:
            return True
        else:
            raise Exception('param_dict has a non-allowed parse_type: %s. Allowed: %s'
                            % (param_dict['parse_type'], allowed_edge_parse_types))
    else:
        raise Exception('required parameter(s) missing : %s'
                        % (required_edge_params - param_dict.keys()))


def _neo4j_sanitize(string):
    """
    Auxilary functions making sure that backdashes are properly noted and managed in Python.

    :param string: string to sanitize
    :return:
    """
    if isinstance(string, str):
        return string.replace('\'', '\"').replace('\\', '\\\\')
    else:
        return string


class GraphDBPipe(object):
    """
    A class that encapsulates the methods needed to work with a database storing the graph of
    biological entities and annotation relationships. This implements a neo4j backend, but can be
    subclassed to use a different backend by overriding _driver static method and the _XXX static
    methods that actually implement what functions XXX() expose.
    """

    def __init__(self):  # REFACTOR: add an option to supply uri, user and pwd not from configs
        password = os.environ['NEOPASS']
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
        if neo4j_db_name is None:
            self._active_database = DEFAULT_DATABASE  # INTEST:
        else:
            self._active_database = neo4j_db_name

    def close(self):
        """
        Closes the connection to the database

        :return:
        """
        try:
            self._driver.close()
        except TypeError:
            pass

    def create(self,
               node_type: db_n_type,
               param_dict: dict) -> Node:
        """
        Creates a node

        :param node_type: neo4j type of the node
        :param param_dict: parameters of the node to be set in the database
        :return: the created node
        """
        with self._driver.session(database=self._active_database) as session:
            new_node = session.write_transaction(self._create, node_type, param_dict)
            return new_node

    @staticmethod
    def _create(tx, node_type, param_dict):
        _check_node_params(param_dict)

        if node_type == 'GOTerm' and param_dict['parse_type'] != 'annotation':
            raise Exception('Term class name and type inconsistency detected')

        instruction_puck = ["CREATE (n:%s)" % node_type]
        set_puck = []

        for key, value in param_dict.items():
            set_puck.append("SET n.%s = '%s'" % (key, _neo4j_sanitize(value)))

        instruction_puck += set_puck
        instruction_puck.append("RETURN n")
        instruction = ' '.join(instruction_puck)
        result = tx.run(instruction)

        return result.single()['n']

    def delete(self,
               node_id: db_id,
               node_type: db_n_type = None) -> Node:
        """
        Deletes a node based on its id and potentially the node type

        :param node_id: internal database id of the node to be deleted
        :param node_type: (optional type of the node to be deleted)
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            # REFACTOR: [better environment]: organism-specific: session(database=db_org)
            deleted = session.write_transaction(self._delete, node_id, node_type)
            return deleted

    @staticmethod
    def _delete(tx, node_id, node_type):
        if node_type is None:
            result = tx.run("MATCH (n) "
                            "WHERE ID(n) = %s "
                            "OPTIONAL MATCH (n)<-[r:annotates]-(a) "
                            "DETACH DELETE a, n " % node_id)
        else:
            result = tx.run("MATCH (n:%s)<-[r:annotates]-(a) "
                            "WHERE ID(n) = %s "
                            "OPTIONAL MATCH (n)<-[r:annotates]-(a) "
                            "DETACH DELETE a, n " % node_type, node_id)
        return [res for res in result]

    def delete_all(self, node_type: db_n_type) -> List[Node]:
        """
        Batch-deletes all they nodes of a given type

        :param node_type: types of the nodes to be deleted
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            suppression = session.write_transaction(self._delete_all, node_type)
            return suppression

    @staticmethod
    def _delete_all(tx, nodetype):
        result = tx.run("MATCH (n:%s) "
                        "OPTIONAL MATCH (n)<-[r:annotates]-(a) "
                        "DETACH DELETE a, n " % nodetype)
        return result

    def clear_database(self) -> None:
        """
        Deletes all the nodes and links between them in the database

        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._clear_database)

    @staticmethod
    def _clear_database(tx):
        node_types = tx.run("MATCH (N) RETURN DISTINCT LABELS(N)")
        node_types = [_type for _type in node_types]
        node_types = [_type["LABELS(N)"][0] for _type in node_types]

        log.info('Clearing the database')
        for node_type in node_types:
            tx.run("MATCH (n:%s) "
                   "OPTIONAL MATCH (n)<-[r:annotates]-(a) "
                   "DETACH DELETE a, n " % node_type)
            log.info('deleted %s' % node_type)

    def get(self,
            node_id: db_id,
            node_type: db_n_type = None) -> Node:
        """
        Gets a single node based on its internal database id and optional type

        :param node_id: internal database id of the node to be retrieved
        :param node_type: (optional) type of the node to be retrieved
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            node = session.write_transaction(self._get, node_type, node_id)
            return node

    @staticmethod
    def _get(tx, node_type, node_id):
        if node_type is None:
            result = tx.run("MATCH (n) "
                            "WHERE ID(n) = %s "
                            "RETURN n" % node_id)

        else:
            result = tx.run("MATCH (n:%s) "
                            "WHERE ID(n) = %s "
                            "RETURN n" % (node_type, node_id))

        node = result.single()

        if node is None:
            return None

        else:
            return node[0]

    def get_all(self, node_type: db_n_type) -> List[Node]:
        """
        Batch-gets all the nodes of a given type

        :param node_type: type of the nodes to get
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            nodes = session.write_transaction(self._get_all, node_type)
            return nodes

    @staticmethod
    def _get_all(tx, node_type):
        result = tx.run("MATCH (n:%s) "
                        "RETURN n" % node_type)
        return [node['n'] for node in result]

    def count(self, node_type: db_n_type) -> int:
        """
        Counts all the nodes of a given type

        :param node_type: type of nodes to count
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            nodes_count = session.write_transaction(self._count, node_type)
            return nodes_count

    @staticmethod
    def _count(tx, node_type):
        result = tx.run("MATCH (n:%s) "
                        "RETURN count(distinct n)" % node_type)

        return result.single()['count(distinct n)']

    def find(self,
             filter_dict: dict,
             node_type: db_n_type = None) -> List[Node]:
        """
        Finds all nodes matching the properties supplied by the filter dictionnary

        :param filter_dict: dictionary containing the parameters and values on which to match
        :param node_type: (optional) type of nodes among which to search
        :return: list of found nodes
        """
        with self._driver.session(database=self._active_database) as session:
            nodes = session.write_transaction(self._find, node_type, filter_dict)
            return nodes

    @staticmethod
    def _find(tx, node_type, filter_dict):
        if node_type is not None:
            instruction_puck = ["MATCH (a:%s)" % node_type]
        else:
            instruction_puck = ["MATCH (a)"]

        where_puck = []

        for key, value in filter_dict.items():
            where_puck.append("a.%s = '%s'" % (key, _neo4j_sanitize(value)))

        where_clause = "WHERE " + ' AND '.join(where_puck) + ' '
        instruction_puck.append(where_clause)
        instruction_puck.append('RETURN a')
        instruction = ' '.join(instruction_puck)
        nodes = tx.run(instruction)

        return [node['a'] for node in nodes]

    def link(self,
             node_id_from: db_id,
             node_id_to: db_id,
             link_type: db_e_type = None,
             params: dict = None) -> List[Relationship]:
        """
        Links two nodes based on their IDs with a link of provided type with provided parameters

        :param node_id_from: internal database id of the node from which the links starts
        :param node_id_to: internal database id of the node to which the links end
        :param link_type: link type
        :param params: provided link parameters
        :return: list containing the the created link object
        """
        with self._driver.session(database=self._active_database) as session:
            link = session.write_transaction(self._link_create, node_id_from, node_id_to, link_type, params)
            return link

    @staticmethod
    def _link_create(tx, node_from, node_to, link_type, params):
        _check_edge_params(params)

        if link_type is None:
            link_type = 'default'

        instructions_puck = ["MATCH (a), (b)",
                             "WHERE ID(a) = %s AND ID(b) = %s" % (node_from, node_to),
                             "CREATE (a)-[r:%s]->(b)" % link_type]

        if params is not None:
            set_puck = []
            for key, value in params.items():
                set_puck.append("SET r.%s = '%s'" % (key, _neo4j_sanitize(value)))
            instructions_puck += set_puck

        instructions_puck.append("RETURN r")
        instruction = ' '.join(instructions_puck)
        rels = tx.run(instruction)

        return [rel['r'] for rel in rels]

    def get_linked(self,
                   node_id: db_id,
                   orientation: e_orientation = 'both',
                   link_type: db_e_type = None,
                   link_param_filter: dict = None) -> List[Node]:
        """
        Get all nodes linked to a given node by a link of a given type

        :param node_id: internal database id of the node from which to start
        :param orientation: 'both'|'in'|'out'. default is 'both'
        :param link_type: (optional) type of links to follow
        :param link_param_filter: (optional) only links whose parameters match this filter will
        be followed
        :return: list of nodes that were found
        """
        with self._driver.session(database=self._active_database) as session:
            results = session.write_transaction(self._get_linked, node_id, orientation, link_type, link_param_filter)
            return results

    @staticmethod
    def _get_linked(tx, node_id, orientation, link_type, link_param_filter):
        instructions_puck = []

        if link_type is None:
            if orientation == 'both':
                instructions_puck.append("MATCH (a)-[r]-(b)")
            if orientation == 'in':
                instructions_puck.append("MATCH (a)<-[r]-(b)")
            if orientation == 'out':
                instructions_puck.append("MATCH (a)-[r]->(b)")
        else:
            if orientation == 'both':
                instructions_puck.append("MATCH (a)-[r:%s]-(b)" % link_type)
            if orientation == 'in':
                instructions_puck.append("MATCH (a)<-[r:%s]-(b)" % link_type)
            if orientation == 'out':
                instructions_puck.append("MATCH (a)-[r:%s]->(b)" % link_type)

        instructions_puck.append("WHERE ID(a) = %s" % node_id)

        if link_param_filter is not None:
            where_puck = []
            for key, value in link_param_filter.items():
                where_puck.append("AND r.%s = '%s'" % (key, _neo4j_sanitize(value)))
            instructions_puck += where_puck

        instructions_puck.append("RETURN b")
        instruction = ' '.join(instructions_puck)
        linked_nodes = tx.run(instruction)

        return [node['b'] for node in linked_nodes]

    def set_attributes(self,
                       node_id: db_id,
                       attributes_dict: dict) -> Node:
        """
        For the node designated by the internal database identifier, sets the properties provided
        by the attributes dictionary

        :param node_id: internal database identifier
        :param attributes_dict: dictionary of properties to set
        :return: edited node
        """
        with self._driver.session(database=self._active_database) as session:
            edited_node = session.write_transaction(self._set_attributes, node_id, attributes_dict)
            return edited_node

    @staticmethod
    def _set_attributes(tx, node_id, attributes_dict):
        instructions_puck = ["MATCH (n) WHERE ID(n) = %s" % node_id]
        for key, value in attributes_dict.items():
            instructions_puck.append("SET n.%s = '%s'" % (key, _neo4j_sanitize(value)))
        instructions_puck.append("RETURN n")
        instruction = ' '.join(instructions_puck)
        log.debug(instruction)
        result = tx.run(instruction)

        return result.single()

    def attach_annotation_tag(self,
                              node_id: db_id,
                              annotation_tag: str,
                              tag_type: db_n_type = None,
                              preferential: bool = False,
                              source: str = 'NA') -> Node:
        """
        For a given node, set a tag with its identifier in a database of a given type.

        :param node_id: id of the node we will be annotating
        :param annotation_tag: the external identifier with witch the node will be annotated
        :param tag_type: (optional) type of the annotation
        :param preferential: (optional) if this is the preferential external identifier
        :param: source: (optional) where the link comes from
        :return: new node containing the annotation
        """
        with self._driver.session(database=self._active_database) as session:
            annot_tag = session.write_transaction(self._attach_annotation_tag,
                                                  node_id, annotation_tag, tag_type,
                                                  preferential, source)
            return annot_tag

    @staticmethod
    def _attach_annotation_tag(tx, node_id, annotation_tag, tag_type, preferential, link_source):
        if tag_type is None:
            tag_type = 'undefined'

        if preferential:
            result = tx.run("MATCH (a) "
                            "WHERE ID(a) = %s "
                            "CREATE (b:Annotation) "
                            "SET b.tag = '%s' "
                            "SET b.type = '%s' "
                            "SET b.parse_type = 'xref' "
                            "SET b.source = '%s' "
                            "CREATE (a)<-[r:annotates]-(b) "
                            "SET r.preferential = True "
                            "SET r.parse_type = 'xref' "
                            "SET r.source = '%s' "
                            "RETURN b" % (node_id,
                                          _neo4j_sanitize(annotation_tag),
                                          _neo4j_sanitize(tag_type),
                                          _neo4j_sanitize(link_source),
                                          _neo4j_sanitize(link_source)))

        else:
            result = tx.run("MATCH (a) "
                            "WHERE ID(a) = %s "
                            "CREATE (b:Annotation) "
                            "SET b.tag = '%s' "
                            "SET b.type = '%s' "
                            "SET b.parse_type = 'xref' "
                            "SET b.source = '%s' "
                            "CREATE (a)<-[r:annotates]-(b) "
                            "SET r.preferential = False "
                            "SET r.parse_type = 'xref' "
                            "SET r.source = '%s' "
                            "RETURN b" % (node_id,
                                          _neo4j_sanitize(annotation_tag),
                                          _neo4j_sanitize(tag_type),
                                          _neo4j_sanitize(link_source),
                                          _neo4j_sanitize(link_source)))

        return result.single()

    def get_from_annotation_tag(self,
                                annotation_tag: str,
                                tag_type: db_n_type = None) -> List[Node]:
        """
        Recovers all node annotated by a given external identifier of a given type

        :param annotation_tag: external identifier
        :param tag_type: (optional) the type of the external identifier
        :return: list of nodes that are annotated by the external identifier
        """
        with self._driver.session(database=self._active_database) as session:
            annotated_nodes = session.write_transaction(self._get_from_annotation_tag, annotation_tag, tag_type)
            return annotated_nodes

    @staticmethod
    def _get_from_annotation_tag(tx, annotation_tag, tag_type):
        annotation_tag = annotation_tag.upper()
        if tag_type is None or tag_type == '':
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' "
                            "RETURN target" % _neo4j_sanitize(annotation_tag))

        else:
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' AND annotnode.type = '%s' "
                            "RETURN target" % (_neo4j_sanitize(annotation_tag),
                                               _neo4j_sanitize(tag_type)))

        pre_return_puck = list(set([node['target'] for node in result]))

        # The issue here is that sometimes the same tag annotates different nodes and we need a way
        # to distinguishing them. Thus we check if we have more than one Uniprot returned. If
        # yes, we re-run the query checking for it

        node_types = [list(node.labels)[0] for node in pre_return_puck]

        if node_types.count('UNIPROT') > 1:
            # we need to look for preferential links:
            if tag_type is None or tag_type == '':
                result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                                "WHERE annotnode.tag = '%s' AND r.preferential = True "
                                "RETURN target" % _neo4j_sanitize(annotation_tag))

            else:
                result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                                "WHERE annotnode.tag = '%s' AND annotnode.type = '%s' AND r.preferential = True "
                                "RETURN target" % (
                                    _neo4j_sanitize(annotation_tag), _neo4j_sanitize(tag_type)))

            return_puck = list(set([node['target'] for node in result]))
            node_types = [list(node.labels)[0] for node in return_puck]

            if node_types.count('UNIPROT') != 1:
                log.debug('Preferential matching failed: for %s, \n \t %s \n \t %s' % (annotation_tag,
                                                                                      return_puck,
                                                                                      pre_return_puck))

            if node_types.count('UNIPROT') == 0:
                return_puck = pre_return_puck

        else:
            return_puck = pre_return_puck

        return return_puck

    def attach_all_node_annotations(self,
                                    node_id: db_id,
                                    annot_type_2_annot_list: dict,
                                    preferential: bool = False,
                                    source: str = 'NA') -> List[Node]:
        # TODO: preferential is supposed to be a list here and needs to be zipped in
        #  enumerate [unused for now]
        """
        Attaches all the external identifiers of given types to a node based on its ID

        :param node_id: internal database id of the node to annotate
        :param annot_type_2_annot_list: {annotation_type: [external identifiers]}
        :param preferential: (optional) if the annotation is preferential
        :param source: (optional: source of the annotation)
        :return: list of all newly created annotation nodes
        """
        with self._driver.session(database=self._active_database) as session:
            tx = session.begin_transaction()
            annot_nodes = []
            for annot_type, annot_list in annot_type_2_annot_list.items():
                if isinstance(annot_list, str):
                    annot_list = [annot_list]
                for annot_tag in annot_list:
                    annot_node = self._attach_annotation_tag(tx, node_id, annot_tag,
                                                             annot_type, preferential, source)
                    annot_nodes.append(annot_node)
            tx.commit()
            return annot_nodes

    def batch_insert(self,
                     type_list: List[db_n_type],
                     param_dicts_list: List[dict],
                     batch_size: int = 1000) -> List[Node]:
        """
        Performs a batch insertions of nodes whose types and parameters are specified by lists

        :param type_list: types of nodes to be inserted into the database
        :param param_dicts_list: parameters of nodes to be inserted into teh database
        :param batch_size: (optional) how many nodes insert per batch
        :return: list of created nodes
        """
        with self._driver.session(database=self._active_database) as session:
            tx = session.begin_transaction()
            new_nodes = []
            for i, (n_type, n_params) in enumerate(zip(type_list, param_dicts_list)):
                new_nodes.append(self._create(tx, n_type, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return new_nodes

    def batch_link(self,
                   id_pairs_list: List[Tuple[db_id, db_id]],
                   type_list: List[db_e_type],
                   param_dicts_list: List[dict],
                   batch_size: int = 1000) -> List[List[Relationship]]:
        """
        Performs a batch link of nodes in the list by the links fo type in the list and assign
        them parameters from the dict

        :param id_pairs_list: pairs nodes to link
        :param type_list: types of the links
        :param param_dicts_list: list of parameters to be assinged to the links
        :param batch_size: (optional) links created in each transaction
        :return: list of set links
        """
        with self._driver.session(database=self._active_database) as session:
            tx = session.begin_transaction()
            new_links = []
            for i, (_from, _to), n_type, n_params in enumerate(zip(id_pairs_list, type_list, param_dicts_list)):
                new_links.append(self._link_create(tx, _from, _to, n_type, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return new_links

    def batch_set_attributes(self,
                             id_list: List[db_id],
                             param_dicts_list: List[dict],
                             batch_size: int = 1000) -> List[Node]:
        """
        Batch sets the attributes of nodes

        :param id_list: list of internal db ids of nodes to set the attributes
        :param param_dicts_list: list of dicts of node attributes to be set
        :param batch_size: (optional) nodes to set parameters in each batch
        :return: list of updated nodes
        """
        with self._driver.session(database=self._active_database) as session:
            tx = session.begin_transaction()
            edited_nodes = []
            for i, (n_id, n_params) in enumerate(zip(id_list, param_dicts_list)):
                edited_nodes.append(self._set_attributes(tx, n_id, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return edited_nodes

    # due to multiple matching of the annotations (aka x-refs, we need to return lists of lists)
    def batch_retrieve_from_annotation_tags(self,
                                            annotation_tags_list: List[str],
                                            annotations_types: List[db_n_type],
                                            batch_size=1000) -> List[List[Node]]:
        """
        Batch retrieve all the nodes annotated by a given list of tags. Returns the list if the

        :param annotation_tags_list: list of external db identifiers
        :param annotations_types: list of types of external db identifiers
        :param batch_size: (optional) nodes to find per batch
        :return: list of lists of found physical entity nodes
        """
        log.info('Batch retrieval started with %s elements' % len(annotation_tags_list))
        with self._driver.session(database=self._active_database) as session:
            tx = session.begin_transaction()
            annotated_nodes = []

            if annotations_types is None or isinstance(annotations_types, str):
                for i, annotation_tag in enumerate(annotation_tags_list):
                    annotated_nodes.append(self._get_from_annotation_tag(tx, annotation_tag,
                                                                  annotations_types))
                    if i % batch_size == 0:
                        log.info('\t %.2f %%' % (float(i) / float(len(annotation_tags_list))*100))
                        tx.commit()
                        tx = session.begin_transaction()

            else:  # we assume it's a list
                for i, (annot_tag, annot_type) in enumerate(zip(annotation_tags_list,
                                                                annotations_types)):
                    annotated_nodes.append(self._get_from_annotation_tag(tx, annot_tag, annot_type))
                    if i % batch_size == 0:
                        log.info('\t %.2f %%' % (float(i) / float(len(annotation_tags_list)) * 100))
                        tx.commit()
                        tx = session.begin_transaction()
            log.info('\t 100 %%')
            tx.commit()
            return annotated_nodes

    def build_indexes(self) -> None:
        """
        Build indexes on the Uniprots and annotations for an accelerated search if they haven't
        been build yet

        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._build_indexes)

    @staticmethod
    def _build_indexes(tx):
        # Because neo4j breaks retro compatibility more often than I rebuild the database
        tx.run("CREATE INDEX IF NOT EXISTS FOR (n:UNIPROT) ON (n.legacyID)")
        tx.run("CREATE INDEX IF NOT EXISTS FOR (n:Annotation) ON (n.tag)")
        tx.run("CREATE INDEX IF NOT EXISTS FOR (n:Annotation) ON (n.tag, n.type)")

    def cross_link_on_xrefs(self,
                            annotation_type: List[db_n_type]) -> List[Node]:
        """
        Connects the nodes that have the same identifiers of the type `annotation_type` with an
        `is_likely_same` bidirectional edge.

        :param annotation_type: source of annotation on which the cross-linking is done
        :return: nodes that has been cross-linked
        """
        with self._driver.session(database=self._active_database) as session:
            dirty_nodes = session.write_transaction(self._cross_link_on_xrefs, annotation_type)
            return dirty_nodes

    @staticmethod
    def _cross_link_on_xrefs(tx, annotation_type):
        result = tx.run("MATCH (a:Annotation)--(n:UNIPROT) "
                        "MATCH (b:Annotation)--(m:UNIPROT) "
                        "WHERE a.tag = b.tag AND a.type = '%s' and m.legacyID <> n.legacyID "
                        "CREATE (m)-[r:is_likely_same]->(n) "
                        "CREATE (m)<-[k:is_likely_same]-(n) "
                        "SET r.linked_on = '%s' "
                        "SET k.linked_on = '%s' " % (annotation_type,
                                                     annotation_type,
                                                     annotation_type))
        return [node for node in result]

    def count_go_annotation_cover(self) -> None:
        """
        Computes and writes the number of UNIOPROT nodes each GO term annotates, directly and
        indirectly, computes the informativity of each GO term and finally the total annotation
        information available on the UNIPROT

        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._count_direct_coverage)
            session.write_transaction(self._count_indirect_coverage)
            session.write_transaction(self._count_up_inf_content)

    @staticmethod
    def _count_direct_coverage(tx):
        tx.run("MATCH (n:UNIPROT)--(a:GOTerm) "
               "WITH a, count(distinct n) as dir_links "
               "SET a.direct_links = dir_links")

    @staticmethod
    def _count_indirect_coverage(tx):
        total_up = tx.run("MATCH (n:UNIPROT) RETURN count(distinct n) as tot_links")
        total_up = total_up.single()['tot_links']

        tx.run("MATCH (n:UNIPROT)-[:is_go_annotation]-(b:GOTerm) "
               "OPTIONAL MATCH (n:UNIPROT)-[:is_go_annotation]->(a:GOTerm)-[:is_a_go*]->(b:GOTerm) "
               "WITH b, count(distinct n) as tot_links "
               "SET b.total_links = tot_links "
               "SET b.information_content = log(toFloat(%s)/toFloat(tot_links))" % (total_up))


    @staticmethod
    def _count_up_inf_content(tx):
        tx.run("MATCH (n:UNIPROT)-[:is_go_annotation]-(b:GOTerm) "
               "OPTIONAL MATCH (n:UNIPROT)-[:is_go_annotation]->(a:GOTerm)-[:is_a_go*]->(b:GOTerm) "
               "WITH n, sum(b.information_content) as tot_inf "
               "SET n.total_information = tot_inf")

    def get_preferential_gene_names(self) -> dict:
        """
        Get a dict that maps the legacy ids of all UNIPROT nodes to their single preferential
        external database identifier

        :return: {UNIPROT.LegacyId: preferential annotation tag}
        """
        with self._driver.session(database=self._active_database) as session:
            name_maps = session.write_transaction(self._get_preferential_gene_names)
            return name_maps

    @staticmethod
    def _get_preferential_gene_names(tx):
        name_maps = tx.run("MATCH (n:UNIPROT)-[r:annotates]-(a:Annotation) "
                           "WHERE r.preferential = True AND a.type = 'UNIPROT_GeneName' "
                           "RETURN n, a")

        name_maps = dict((res['n']['legacyID'], res['a']['tag'])
                         for res in name_maps)

        return name_maps

    def check_connection_permutation(self,
                                     legacy_id_1: str,
                                     legacy_id_2: str) -> int:
        """
        Counts the number of paths between Uniprots with given legacy IDs that includes an
        `is_likely_same` edge connection.

        :param legacy_id_1: LegacyId of a first Uniprot node
        :param legacy_id_2: LegacyId of a second Uniprot node
        :return: how many paths including an `is_likely_same` edge are between the two nodes
        """
        with self._driver.session(database=self._active_database) as session:
            legal = session.write_transaction(self._check_connection_permutation, legacy_id_1, legacy_id_2)
            log.info('checking_permutation')
            return legal

    @staticmethod
    def _check_connection_permutation(tx, legacy_id_1, legacy_id_2):
        name_maps_1 = tx.run("MATCH (n:UNIPROT)--(k:UNIPROT)-[:is_likely_same]-(b) "
                            "WHERE (n.legacyID = '%s' AND b.legacyID = '%s') "
                            "OR (n.legacyID = '%s' AND b.legacyID = '%s') "
                            "RETURN n, k, b" % (legacy_id_1, legacy_id_2, legacy_id_2, legacy_id_1))

        legal = len([result for result in name_maps_1])

        if not legal:
            name_maps_2 = tx.run("MATCH (n:UNIPROT)-[:is_likely_same]-(u:UNIPROT)--(k:UNIPROT)-[:is_likely_same]-(b) "
                                "WHERE (n.legacyID = '%s' AND b.legacyID = '%s') "
                                "OR (n.legacyID = '%s' AND b.legacyID = '%s') "
                                "RETURN n, u, k, b" % (legacy_id_1, legacy_id_2, legacy_id_2, legacy_id_1))
            legal = len([result for result in name_maps_2])

        return legal

    def node_stats(self) -> None:
        """
        Prints to log.info the number of nodes in each type present in the database

        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._node_stats)

    @staticmethod
    def _node_stats(tx):
        node_types = tx.run("MATCH (N) RETURN DISTINCT LABELS(N)")
        node_types = [_type for _type in node_types]
        node_types = [_type["LABELS(N)"][0] for _type in node_types]

        for node_type in node_types:
            total = tx.run("MATCH (N:%s) RETURN COUNT(N)" % node_type).single()["COUNT(N)"]
            log.info("%s : %d" % (node_type, total))

    def mark_forbidden_nodes(self, excluded_names_or_leg_ids) -> List[Node]:
        """
        Marks the nodes that are not allowed to be used for information flow computation (too
        high degree or too high eigenvector weight)

        :param excluded_names_or_leg_ids: list of legacy ID or names of nodes to exclude
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            nodes_list = session.write_transaction(self._mark_forbidden_nodes,
                                                   excluded_names_or_leg_ids)
        return nodes_list

    @staticmethod
    def _mark_forbidden_nodes(tx, excluded_names_or_leg_ids):
        all_nodes = []
        for tag in excluded_names_or_leg_ids:
            nodes = tx.run("MATCH (N) "
                           "WHERE (N.displayName='%s' OR N.legacyID='%s') "
                           "AND N.parse_type='physical_entity' "
                           "SET N.forbidden='True' "
                           "RETURN N" % (tag, tag))
            nodes = [_node['N'] for _node in nodes]
            all_nodes += nodes

        return all_nodes

    def parse_physical_entity_net(self, main_connex_only: bool = False) -> (Dict[int, Node],
                                                            List[Relationship]):
        """
        Runs and prints diagnostics on its own content

        :return:
        """
        log.info('Massive pull from the database, this might take a while, please wait')
        with self._driver.session(database=self._active_database) as session:
            nodes_dict, rels_list = session.write_transaction(self._parse_physical_entity_net,
                                                              main_connex_only)
        log.info('Pull suceeded')
        return nodes_dict, rels_list

    @staticmethod
    def _parse_physical_entity_net(tx, main_connex_only):

        if main_connex_only:
            nodes = tx.run(
                "MATCH (N)-[r]-(M) "
                "WHERE (N.parse_type='physical_entity' AND M.parse_type='physical_entity') "
                "AND N.main_connex='True' "
                "AND NOT (EXISTS(N.forbidden) OR EXISTS(M.forbidden)) "
                "AND (r.parse_type='physical_entity_molecular_interaction' "
                "OR r.parse_type='identity' OR r.parse_type='refines')"
                "WITH [N, M] as tl "
                "UNWIND tl as n "
                "RETURN DISTINCT(n)")

            nodes_dict = {_node['n'].id: _node['n'] for _node in nodes}

            rels = tx.run(
                "MATCH (N)-[r]-(M) "
                "WHERE (N.parse_type='physical_entity' AND M.parse_type='physical_entity') "
                "AND N.main_connex='True' "
                "AND NOT (EXISTS(N.forbidden) OR EXISTS(M.forbidden)) "
                "AND (r.parse_type='physical_entity_molecular_interaction' "
                "OR r.parse_type='identity' OR r.parse_type='refines')"
                "RETURN DISTINCT(r)")

            rels_list = [_rel['r']for _rel in rels]

            return nodes_dict, rels_list

        else:
            nodes = tx.run(
                "MATCH (N)-[r]-(M) "
                "WHERE (N.parse_type='physical_entity' AND M.parse_type='physical_entity') "
                "AND NOT (EXISTS(N.forbidden) OR EXISTS(M.forbidden)) "
                "AND (r.parse_type='physical_entity_molecular_interaction' "
                "OR r.parse_type='identity' OR r.parse_type='refines')"
                "WITH [N, M] as tl "
                "UNWIND tl as n "
                "RETURN DISTINCT(n)")

            nodes_dict = {_node['n'].id: _node['n'] for _node in nodes}

            rels = tx.run(
                "MATCH (N)-[r]-(M) "
                "WHERE (N.parse_type='physical_entity' AND M.parse_type='physical_entity') " 
                "AND NOT (EXISTS(N.forbidden) OR EXISTS(M.forbidden)) "
                "AND (r.parse_type='physical_entity_molecular_interaction' "
                "OR r.parse_type='identity' OR r.parse_type='refines')"
                "RETURN DISTINCT(r)")

            rels_list = [_rel['r']for _rel in rels]

            return nodes_dict, rels_list

    def parse_knowledge_entity_net(self) -> (Dict[int, Node], List[Relationship]):
        """
        Runs and prints diagnostics on its own content

        :return:
        """
        log.info('Massive pull from the database, this might take a while, please wait')
        with self._driver.session(database=self._active_database) as session:
            nodes_dict, rels_list = session.write_transaction(self._parse_knowledge_entity_net)
        log.info('Pull suceeded')
        return nodes_dict, rels_list

    @staticmethod
    def _parse_knowledge_entity_net(tx):

        # This is directional. If matched, uniprot or physical entity will always be first.
        nodes = tx.run(
            "MATCH (N)-[r]-(M) "
            "WHERE ((N.parse_type='physical_entity' OR N.parse_type='annotation') "
            "AND M.parse_type='annotation') "
            "WITH [N, M] as tl "
            "UNWIND tl as n "
            "RETURN DISTINCT(n)")

        nodes_dict = {_node['n'].id: _node['n'] for _node in nodes}

        rels = tx.run(
            "MATCH (N)-[r]->(M) "
            "WHERE ((N.parse_type='physical_entity' OR N.parse_type='annotation') "
            "AND M.parse_type='annotation') "
            "RETURN DISTINCT(r)")

        rels_list = [_rel['r']for _rel in rels]

        return nodes_dict, rels_list


    def erase_node_properties(self, properties_list):
        """
        Removes all proprties whose name is in the list from all the nodes who have it.

        :param properties_list: list of properties to be removed
        :raise Exception: if a parameter that is required is cleared
        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._erase_node_properties, properties_list)

    @staticmethod
    def _erase_node_properties(tx, properties_list):
        for _property in properties_list:
            if _property in required_node_params:
                raise Exception('Trying to clear a required node parameter %s' % _property)
            resets = tx.run("MATCH (N) "
                            "WHERE EXISTS(N.%s) "
                            "REMOVE N.%s "
                            "RETURN COUNT(N)" % (_property, _property)).single()['COUNT(N)']
            log.info("Cleared property %s in %d nodes" % (_property, resets))

    def self_diag(self) -> None:
        """
        Runs and prints diagnostics on its own content to the stdout

        :return:
        """
        with self._driver.session(database=self._active_database) as session:
            session.write_transaction(self._self_diag)

    @staticmethod
    def _self_diag(tx):  # REFACTOR: split. high cyclomatic complexity (~18)
        node_types = tx.run("MATCH (N) RETURN DISTINCT LABELS(N)")
        node_types = [_type for _type in node_types]
        node_types = [_type["LABELS(N)"][0] for _type in node_types]
        # TODO: add multi-label support
        rel_types = tx.run("MATCH ()-[r]-() RETURN DISTINCT TYPE(r)")
        rel_types = [_type["TYPE(r)"] for _type in rel_types]

        print('\nNodes characterization:')

        for node_type in node_types:

            total = tx.run("MATCH (N:%s) RETURN COUNT(N)" % node_type).single()["COUNT(N)"]

            print("%s: %d" % (node_type, total))

            distinct_keysets = tx.run("MATCH (N:%s) RETURN DISTINCT keys(N)" % node_type)
            # distinct_keysets = set([frozenset(keyset) in distinct_keysets])
            distinct_properties_for_type = set()
            distinct_keysets = [keyset["keys(N)"] for keyset in distinct_keysets]
            [distinct_properties_for_type.update(keyset) for keyset in distinct_keysets]

            for _property in distinct_properties_for_type:

                property_distinct_values = tx.run(
                    "MATCH (N:%s) "
                    "RETURN COUNT(DISTINCT N.%s)"
                    % (node_type, _property)).single()["COUNT(DISTINCT N.%s)" % _property]

                property_occurences = tx.run(
                    "MATCH (N:%s) "
                    "WHERE EXISTS(N.%s) "
                    "RETURN COUNT(N)" % (node_type, _property)).single()["COUNT(N)"]

                print("\t %s: %d/%.2f %%; %d distinct values"
                      % (_property, property_occurences,
                         property_occurences/total*100, property_distinct_values))

                if property_distinct_values < 6:

                    property_distinct_values = tx.run("MATCH (N:%s) RETURN DISTINCT N.%s"
                                                      % (node_type, _property))

                    for _value in property_distinct_values:
                        _val = _value[0]
                        _val_occurences  = tx.run(
                            "MATCH (N:%s) "
                            "WHERE N.%s='%s' "
                            "RETURN COUNT(N)" % (node_type, _property, _val)).single()[0]


                        print("\t\t %s: %d/%.2f %%"
                              % (_val, _val_occurences, _val_occurences/property_occurences*100))


        print('\nEdges characterization:')

        for rel_type in rel_types:

            total = tx.run("MATCH ()-[r:%s]-() RETURN COUNT(r)" % rel_type).single()["COUNT(r)"]

            print("%s: %d" % (rel_type, total))

            distinct_keysets = tx.run("MATCH ()-[r:%s]-() RETURN DISTINCT keys(r)" % rel_type)
            distinct_properties_for_type = set()
            distinct_keysets = [keyset["keys(r)"] for keyset in distinct_keysets]
            [distinct_properties_for_type.update(keyset) for keyset in distinct_keysets]

            for _property in distinct_properties_for_type:
                property_distinct_values = tx.run("MATCH ()-[r:%s]-() "
                                                  "RETURN COUNT(DISTINCT r.%s)"
                                                  % (rel_type, _property)).single()[0]

                property_occurences = tx.run("MATCH ()-[r:%s]-() "
                                             "WHERE EXISTS(r.%s) "
                                             "RETURN COUNT(r)" % (rel_type, _property)).single()[0]

                print("\t %s: %d/%.2f %%; %d distinct values"
                      % (_property, property_occurences,
                         property_occurences/total*100, property_distinct_values))

                if property_distinct_values < 6:
                    property_distinct_values = tx.run("MATCH ()-[r:%s]-() "
                                                      "RETURN DISTINCT r.%s"
                                                      % (rel_type, _property))

                    for _value in property_distinct_values:
                        _val = _value[0]
                        _val_occurences = tx.run(
                            "MATCH ()-[r:%s]-() "
                            "WHERE r.%s='%s' "
                            "RETURN COUNT(r)"
                            % (rel_type, _property, _val)).single()[0]


                        print("\t\t %s: %d/%.2f %%"
                              % (_val, _val_occurences, _val_occurences/property_occurences*100))

        print('\nPattern occurences:')
        for A_type, B_type in combinations_with_replacement(node_types, 2):
            for r_type in rel_types:
                occurences = tx.run("MATCH (A:%s)-[r:%s]-(B:%s) "
                                    "RETURN COUNT(r) "
                                    "AS occurences"
                                    % (A_type, r_type, B_type)).single()['occurences']

                if occurences:
                    print("\t(%s)-%s-(%s) : %d" % (A_type, r_type, B_type, occurences))

        print('\nParse_type pattern occurences:')
        for A_parse_type, B_parse_type in combinations_with_replacement(allowed_node_parse_types, 2):
            for r_parse_type in allowed_edge_parse_types:
                occurences = tx.run(
                    "MATCH (A {parse_type:'%s'})-[r {parse_type:'%s'}]-(B {parse_type:'%s'}) "
                    "RETURN COUNT(r) "
                    "AS occurences"
                    % (A_parse_type, r_parse_type, B_parse_type)).single()['occurences']

                if occurences:
                    print("\t(%s)-%s-(%s) : %d"
                          % (A_parse_type, r_parse_type, B_parse_type, occurences))


        print('\nHighest degree nodes:')
        node_2_degree_map = tx.run("MATCH (N {parse_type: 'physical_entity'})-"
                                   "-(M {parse_type: 'physical_entity'})"
                                   "RETURN N,COUNT(M)")

        node_2_degree_map = [(_match['N'], int(_match['COUNT(M)'])) for _match in node_2_degree_map]
        node_2_degree_map = sorted(node_2_degree_map, key=lambda x: x[1])[100:]
        node_2_degree_map = ['%s, %s: %d' % (node['legacyID'], node.id, count)
                             for node, count in node_2_degree_map]

        for entry in node_2_degree_map:
            print(entry)



if __name__ == "__main__":
    neo4j_pipe = GraphDBPipe()
    # nodes_dict, edges_dict = neo4j_pipe.parse_physical_entity_net()
    # print(len(nodes_dict.keys()), len(edges_dict.keys()))
    # print(next(iter(nodes_dict.items())), next(iter(edges_dict.items())))
    neo4j_pipe.self_diag()

    # neo4j_pipe.delete_all('Test')
    # neo4j_pipe.delete_all('RNA')
    # neo4j_pipe.delete_all('Annotation')
    # neo4j_pipe.delete_all('Location')
    # neo4j_pipe.delete_all('Protein')
    # print(neo4j_pipe.create('Protein', {"chat": "miau", "dog": "waf"}))
    # print neo4j_pipe.create('Complex', {"chat": "miau", "dog": "waf"})
    # neo4j_pipe.link(1, 3, 'test', {"weight": 2, "source": "Andrei"})
    # print neo4j_pipe.get(0)['custom']
    # print neo4j_pipe.get_all("Protein")[:10]
    # print neo4j_pipe.delete_all('Complex')
    # print neo4j_pipe.count("Annotation")
    # neo4j_pipe.get_linked(1, 'both', 'test', {"source": "Andrei"})
    # neo4j_pipe.delete(3, 'Protein')
    # print neo4j_pipe.get(1, "Complex").id
    # result_list = neo4j_pipe.find({"chat": "miau"})
    # print result_list
    # print [node.id for node in result_list]
    # neo4j_pipe.set_attributes(1, {"bird": "cookoo"})
    # print neo4j_pipe.attach_annotation_tag(44, "Q123541", "UP Acc ID")
    # nodes = neo4j_pipe.get_from_annotation_tag("MPP4", None)
    # print nodes
    # super_nodes = neo4j_pipe.batch_retrieve_from_annotation_tags(['RNF14', 'REM1', 'GAG', 'MZT2A', 'FHIT'], None)
    # for node_set in super_nodes:
    #     log.info(node_set)
    # print nodes[0]._values[0]['dog']
    # print neo4j_pipe.attach_all_node_annotations(0, {'super': ['wo1', 'wo2'], 'super-duper': ['aka', 'coco']})
    # neo4j_pipe.batch_insert(["Protein", "Complex"], [{'a': 1}, {'a': 2, 'b': 3}])
    # neo4j_pipe.batch_link([(25, 26), (25, 1)], [None, 'reaction'], [None, {"weight": 3, "source": "test"}])
    # print neo4j_pipe.batch_set_attributes([0, 1, 2, 44, 45], [{'custom': 'Main_connex'}]*5)

    # print neo4j_pipe.cross_link_on_xrefs()

    # print neo4j_pipe.count_go_annotation_cover()

