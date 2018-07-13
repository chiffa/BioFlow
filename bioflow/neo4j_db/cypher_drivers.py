import os
from pprint import pprint
from neo4j.v1 import GraphDatabase
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)

uri = 'bolt://localhost:7687'
user = 'neo4j'
password = os.environ['NEOPASS']


# TODO: add the instructions echo
class GraphDBPipe(object):

    def __init__(self):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def create(self, node_type, param_dict):
        with self._driver.session() as session:
            new_node = session.write_transaction(self._create, node_type, param_dict)
            return new_node

    @staticmethod
    def _create(tx, node_type, param_dict):
        instruction_puck = ["CREATE (n:%s)" % node_type]
        set_puck = []

        for key, value in param_dict.iteritems():
            set_puck.append("SET n.%s = '%s'" % (key, value))

        instruction_puck += set_puck
        instruction_puck.append("RETURN n")
        instruction = ' '.join(instruction_puck)
        result = tx.run(instruction)

        return result.single()

    def delete(self, node_id, node_type=None):
        with self._driver.session() as session:
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

    def delete_all(self, node_type):
        with self._driver.session() as session:
            supression = session.write_transaction(self._delete_all, node_type)
            return supression

    @staticmethod
    def _delete_all(tx, nodetype):
        result = tx.run("MATCH (n:%s) "
                        "OPTIONAL MATCH (n:)<-[r:annotates]-(a) "
                        "DETACH DELETE a, n " % nodetype)
        return result

    def get(self, node_id, node_type=None):
        with self._driver.session() as session:
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

    def get_all(self, node_type):
        with self._driver.session() as session:
            nodes = session.write_transaction(self._get_all, node_type)
            return nodes

    @staticmethod
    def _get_all(tx, node_type):
        result = tx.run("MATCH (n:%s) "
                        "RETURN n" % node_type)
        return [node['n'] for node in result]

    def count(self, node_type):
        with self._driver.session() as session:
            nodes_count = session.write_transaction(self._count, node_type)
            return nodes_count

    @staticmethod
    def _count(tx, node_type):
        result = tx.run("MATCH (n:%s) "
                        "RETURN count(distinct n)" % node_type)

        return result.single()['count(distinct n)']

    def find(self, filter_dict, node_type=None):
        with self._driver.session() as session:
            nodes = session.write_transaction(self._find, node_type, filter_dict)
            return nodes

    @staticmethod
    def _find(tx, node_type, filter_dict):
        if node_type is not None:
            instruction_puck = ["MATCH (a:%s)" % node_type]
        else:
            instruction_puck = ["MATCH (a)"]

        where_puck = []

        for key, value in filter_dict.iteritems():
            where_puck.append("a.%s = '%s'" % (key, value))

        where_clause = "WHERE " + ' AND '.join(where_puck) + ' '
        instruction_puck.append(where_clause)
        instruction_puck.append('RETURN a')
        instruction = ' '.join(instruction_puck)
        nodes = tx.run(instruction)

        return [node['a'] for node in nodes]

    def link(self, node_id_from, node_id_to, link_type=None, params=None):
        with self._driver.session() as session:
            link = session.write_transaction(self._link_create, node_id_from, node_id_to, link_type, params)
            return link

    @staticmethod
    def _link_create(tx, node_from, node_to, link_type, params):

        if link_type is None:
            link_type = 'default'

        instructions_puck = ["MATCH (a), (b)",
                             "WHERE ID(a) = %s AND ID(b) = %s" % (node_from, node_to),
                             "CREATE (a)-[r:%s]->(b)" % link_type]

        if params is not None:
            set_puck = []
            for key, value in params.iteritems():
                set_puck.append("SET r.%s = '%s'" % (key, value))
            instructions_puck += set_puck

        instructions_puck.append("RETURN r")
        instruction = ' '.join(instructions_puck)
        rels = tx.run(instruction)

        return [rel['r'] for rel in rels]

    def get_linked(self, node_id, orientation='both', link_type=None, link_param_filter=None):
        with self._driver.session() as session:
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
            for key, value in link_param_filter.iteritems():
                where_puck.append("AND r.%s = '%s'" % (key, value))
            instructions_puck += where_puck

        instructions_puck.append("RETURN b")
        instruction = ' '.join(instructions_puck)
        linked_nodes = tx.run(instruction)

        return [node['b'] for node in linked_nodes]

    def set_attributes(self, node_id, attributes_dict):
        with self._driver.session() as session:
            edited_node = session.write_transaction(self._set_attributes, node_id, attributes_dict)
            return edited_node

    @staticmethod
    def _set_attributes(tx, node_id, attributes_dict):
        instructions_puck = ["MATCH (n) WHERE ID(n) = %s" % node_id]
        for key, value in attributes_dict.iteritems():
            instructions_puck.append("SET n.%s = '%s'" % (key, value))
        instructions_puck.append("RETURN n")
        instruction = ' '.join(instructions_puck)
        log.debug(instruction)
        result = tx.run(instruction)

        return result.single()

    def attach_annotation_tag(self, node_id, annotation_tag, tag_type=None):
        with self._driver.session() as session:
            annot_tag = session.write_transaction(self._attach_annotation_tag, node_id, annotation_tag, tag_type)
            return annot_tag

    @staticmethod
    def _attach_annotation_tag(tx, node_id, annotation_tag, tag_type):
        if tag_type is None:
            tag_type = 'undefined'

        result = tx.run("MATCH (a) "
                        "WHERE ID(a) = %s "
                        "CREATE (b:Annotation) "
                        "SET b.tag = '%s' "
                        "SET b.type = '%s' "
                        "CREATE (a)<-[r:annotates]-(b) "
                        "RETURN b" % (node_id, annotation_tag, tag_type))

        return result.single()

    def get_from_annotation_tag(self, annotation_tag, tag_type=None):
        with self._driver.session() as session:
            annotated_nodes = session.write_transaction(self._get_from_annotation_tag, annotation_tag, tag_type)
            return annotated_nodes

    @staticmethod
    def _get_from_annotation_tag(tx, annotation_tag, tag_type):
        if tag_type is None:
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' "
                            "RETURN target" % annotation_tag)

        else:
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' AND annotnode.type = '%s' "
                            "RETURN target" % (annotation_tag, tag_type))

        return list(set([node['target'] for node in result]))

    def attach_all_node_annotations(self, node_id, annot_type_2_annot_list):
        with self._driver.session() as session:
            tx = session.begin_transaction()
            annot_nodes = []
            for annot_type, annot_list in annot_type_2_annot_list.iteritems():
                for annot_tag in annot_list:
                    annot_node = self._attach_annotation_tag(tx, node_id, annot_tag, annot_type)
                    annot_nodes.append(annot_node)
            tx.commit()
            return annot_nodes

    def batch_insert(self, type_list, param_dicts_list, batch_size=1000):
        with self._driver.session() as session:
            tx = session.begin_transaction()
            new_nodes = []
            for i, (n_type, n_params) in enumerate(zip(type_list, param_dicts_list)):
                new_nodes.append(self._create(tx, n_type, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return new_nodes

    def batch_link(self, id_pairs_list, type_list, param_dicts_list, batch_size=1000):
        with self._driver.session() as session:
            tx = session.begin_transaction()
            new_links = []
            for i, (_from, _to), n_type, n_params in enumerate(zip(id_pairs_list, type_list, param_dicts_list)):
                new_links.append(self._link_create(tx, _from, _to, n_type, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return new_links

    def batch_set_attributes(self, id_list, param_dicts_list, batch_size=1000):
        with self._driver.session() as session:
            tx = session.begin_transaction()
            edited_nodes = []
            for i, (n_id, n_params) in enumerate(zip(id_list, param_dicts_list)):
                edited_nodes.append(self._set_attributes(tx, n_id, n_params))
                if i % batch_size == 0:
                    tx.commit()
                    tx = session.begin_transaction()
            tx.commit()
            return edited_nodes


if __name__ == "__main__":
    neo4j_pipe = GraphDBPipe()
    # neo4j_pipe.delete_all('Test')
    # neo4j_pipe.create('Protein', {"chat": "miau", "dog": "waf"})
    # print neo4j_pipe.create('Complex', {"chat": "miau", "dog": "waf"})
    # neo4j_pipe.link(1, 3, 'test', {"weight": 2, "source": "Andrei"})
    # print neo4j_pipe.get(0).properties['custom']
    # print neo4j_pipe.get_all("Protein")
    # print neo4j_pipe.delete(44) # WTF?
    # print neo4j_pipe.count("Annotation")
    # neo4j_pipe.get_linked(1, 'both', 'test', {"source": "Andrei"})
    # neo4j_pipe.delete(3, 'Protein')
    # print neo4j_pipe.get(1, "Complex").id
    # result_list = neo4j_pipe.find({"chat": "miau"})
    # print result_list
    # print [node.id for node in result_list]
    # neo4j_pipe.set_attributes(1, {"bird": "cookoo"})
    # print neo4j_pipe.attach_annotation_tag(44, "Q123541", "UP Acc ID")
    # nodes = neo4j_pipe.get_from_annotation_tag("Q123541", "UP Acc ID")
    # print nodes[0]._values[0].properties['dog']
    # print neo4j_pipe.attach_all_node_annotations(0, {'super': ['wo1', 'wo2'], 'super-duper': ['aka', 'coco']})
    # neo4j_pipe.batch_insert(["Protein", "Complex"], [{'a': 1}, {'a': 2, 'b': 3}])
    # neo4j_pipe.batch_link([(25, 26), (25, 1)], [None, 'reaction'], [None, {"weight": 3, "source": "test"}])
    # print neo4j_pipe.batch_set_attributes([0, 1, 2, 44, 45], [{'custom': 'Main_connex'}]*5)

