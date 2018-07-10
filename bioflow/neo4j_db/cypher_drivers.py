from neo4j.v1 import GraphDatabase
import code


def interact():
    code.InteractiveConsole(locals=globals()).interact()

uri = 'bolt://localhost:7687'
user = 'neo4j'
password = ''

# the question we need to answer is whether we want to capture the param_dict from the class properties
# derived from that one

# and I am not sure it is a good idea to keep imitating the bulbs neo4j node = python object.
# perhaps a meta-driver will be sufficient?


class GraphDBPipe(object):

    def __init__(self):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def create(self, node_type, param_dict):
        with self._driver.session as session:
            new_node = session.write_transaction(self._create, node_type, param_dict)
            print new_node

    @staticmethod
    def _create(tx, node_type, param_dict):
        instruction_puck = ["CREATE (n:%s)" % node_type]
        set_puck = []

        for key, value in param_dict:
            set_puck.append("SET n.%s = '%s'" % (key, value))

        instruction_puck += set_puck
        instruction_puck.append("RETURN n")
        instruction = ' '.join(instruction_puck)

        result = tx.run(instruction)
        # debugging
        interact()
        return result.single()

    def delete(self, node_id, node_type=None):
        with self._driver.session as session:
            deleted = session.write_transaction(self._delete, node_id, node_type)
            print deleted

    @staticmethod
    def _delete(tx, node_id, node_type):
        if node_type is None:
            result = tx.run("MATCH (n) "
                            "WHERE ID(n) = %s"
                            "DETACH DELETE n " % node_id)
        else:
            result = tx.run("MATCH (n:%s) "
                            "WHERE ID(n) = %s"
                            "DETACH DELETE n " % (node_type, node_id))
        return result

    def delete_all(self, node_type):
        with self._driver.session() as session:
            supression = session.write_transaction(self._delete_all, node_type)
            print supression

    @staticmethod
    def _delete_all(tx, nodetype):
        result = tx.run("MATCH (n:%s) "
                        "DETACH DELETE n " % nodetype)
        return result

    def get(self, node_id):
        with self._driver.session() as session:
            node = session.write_transaction(self._get, self._type, node_id)
            print node

    @staticmethod
    def _get(tx, node_type, node_id):
        result = tx.run("MATCH (a:%s) "
                        "WHERE ID(a) = %s "
                        "RETURN a" % (node_type, node_id))

        return result.single()[0]

    def find(self, filter_dict, node_type=None):
        with self._driver.session() as session:
            nodes = session.write_transaction(self._find, node_type, filter_dict)
            print nodes

    @staticmethod
    def _find(tx, node_type, filter_dict):
        if node_type is not None:
            instruction_puck = ["MATCH (a:%s)" % node_type]
        else:
            instruction_puck = ["MATCH (n)"]

        where_puck = []

        for key, value in filter_dict.iteritems():
            where_puck.append("a.%s = '%s'" % (key, value))

        where_clause = "WHERE " + ' AND '.join(where_puck) + ' '
        instruction_puck.append(where_clause)
        instruction_puck.append('RETURN a')
        instruction = ' '.join(instruction_puck)
        nodes = tx.run(instruction)

        return [node for node in nodes]

    def link(self, node_id_from, node_id_to, link_type=None, params=None):
        with self._driver.session() as session:
            link = session.write_transaction(self._link_create, node_id_from, node_id_to, link_type, params)
            print link

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
        rel = tx.run(instruction)

        return rel

    def get_linked(self, node_id, orientation='both', link_type=None, link_param_filter=None):
        with self._driver.session() as session:
            results = session.write_transaction(self._get_linked, node_id, orientation, link_type, link_param_filter)
            print results

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
                where_puck.append("WHERE r.%s = '%s'" % (key, value))
            instructions_puck += where_puck

        instructions_puck.append("RETURN b")
        instruction = ' '.join(instructions_puck)
        linked_nodes = tx.run(instruction)

        return linked_nodes

    def set_attributes(self, node_id, attributes_dict):
        with self._driver.session() as session:
            edited_node = session.write_transaction(self.cypher_edit_node, node_id, attributes_dict)
            print edited_node

    @staticmethod
    def cypher_edit_node(tx, node_id, attributes_dict):
        instructions_puck = ["MATCH (n)"]
        instructions_puck.append("WHERE ID(n) = %s" % node_id)
        for key, value in attributes_dict:
            instructions_puck.append("SET n.%s = '%s'" % (key, value))
        instructions_puck.append("RETURN n")
        instruction = ' '.join(instructions_puck)
        result = tx.run(instruction)

        return result.single()

    def attach_annotation_tag(self, node_id, annotation_tag, tag_type=None):
        with self._driver.session() as session:
            annot_tag = session.write_transaction(self.cypher_attach_annotation_tag, node_id, annotation_tag, tag_type)
            print annot_tag

    @staticmethod
    def cypher_attach_annotation_tag(tx, node_id, annotation_tag, tag_type):
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
            annotated_nodes = session.write_transaction(self.cypher_get_from_annotation_tag, annotation_tag, tag_type)
            print annotated_nodes

    @staticmethod
    def cypher_get_from_annotation_tag(tx, annotation_tag, tag_type):
        if tag_type is None:
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' "
                            "RETURN target" % annotation_tag)

        else:
            result = tx.run("MATCH (annotnode:Annotation)-[r:annotates]->(target) "
                            "WHERE annotnode.tag = '%s' AND annotnode.type = '%s' "
                            "RETURN target" % (annotation_tag, tag_type))

        return [node for node in result]


class TestNode(object):

    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
        self.node_type = 'Test'

    def close(self):
        self._driver.close()

    def delete_all(self):
        with self._driver.session() as session:
            supression = session.write_transaction(self.cypher_delete_all, self.node_type)
            print supression

    @staticmethod
    def cypher_delete_all(tx, nodetype):
        result = tx.run("MATCH (n:%s) "
                        "DETACH DELETE n " % nodetype)

        return result

    def create_node(self):
        with self._driver.session() as session:
            new_node = session.write_transaction(self.cypher_create_node, self.node_type)
            print new_node

    @staticmethod
    def cypher_create_node(tx, nodetype):
        result = tx.run("CREATE (a:%s) "
                        "SET a.nodetype = '%s' "
                        "RETURN a.nodetype + ' ' + id(a)" % (nodetype, nodetype))
        return result.single()[0]

    def set_attribute(self, node_id, attribute_name, attribute_value):
        with self._driver.session() as session:
            edited_node = session.write_transaction(self.cypher_edit_node, node_id, attribute_name, attribute_value)
            print edited_node

    @staticmethod
    def cypher_edit_node(tx, node_id, attribute_name, attribute_value):
        result = tx.run("MATCH (n) "
                        "WHERE ID(n) = %s "
                        "SET n.%s = '%s' "
                        "RETURN n" % (node_id, attribute_name, attribute_value))
        return result.single()[0]

    def get_linked(self, node_id, link_type=None):
        with self._driver.session() as session:
            linked = session.write_transaction(self.cypher_get_linked, node_id, link_type)
            print linked

    @staticmethod
    def cypher_get_linked(tx, node_id, link_type):
        if link_type is None:
            result = tx.run("MATCH (a)--(b) "
                            "WHERE ID(a) = %s "
                            "RETURN b" % node_id)

        else:
            result = tx.run("MATCH (a)-[r:%s]-(b) "
                            "WHERE ID(a) = %s "
                            "RETURN b" % (link_type, node_id))

        return [node for node in result]

    def link(self, node_id_1, node_id_2, link_type=None):
        with self._driver.session() as session:
            linked = session.write_transaction(self.cypher_link, node_id_1, node_id_2, link_type)
            print linked

    @staticmethod
    def cypher_link(tx, node_id_1, node_id_2, link_type):
        if link_type is None:
            link_type = 'default'

        result = tx.run("MATCH (a), (b) "
                        "WHERE ID(a) = %s AND ID(b) = %s "
                        "CREATE (a)-[r:%s]->(b) "
                        "RETURN r" % (node_id_1, node_id_2, link_type))

        return result.single()




if __name__ == "__main__":
    # neo4j_instance = HelloWorldExample('bolt://localhost:7687', 'neo4j', 'sergvcx')
    # neo4j_instance.print_greetings('Hello World!')

    # neo4j_instance = TestNode('bolt://localhost:7687', 'neo4j', 'sergvcx')
    # neo4j_instance.node_type = 'Greeting'
    # neo4j_instance.create_node()
    # neo4j_instance.create_node()
    # neo4j_instance.create_node()
    # neo4j_instance.link(2, 22)
    # neo4j_instance.link(2, 23, 'test')
    # neo4j_instance.get_linked(2)
    # neo4j_instance.get_linked(2, 'test')
    # neo4j_instance.attach_annotation_tag(2, 'TEST_ANNOT_TAG', 'TEST_TYPE')
    # neo4j_instance.get_from_annotation_tag('TEST_ANNOT_TAG', 'TEST_TYPE')
    # neo4j_instance.delete_all()
    # neo4j_instance.set_attribute(23, 'kitty', 'meow')

    neo4j_instance = GraphDBPipe()
