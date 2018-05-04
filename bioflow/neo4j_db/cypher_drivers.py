from neo4j.v1 import GraphDatabase


class HelloWorldExample(object):

    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def print_greetings(self, message):
        with self._driver.session() as session:
            greeting = session.write_transaction(self._create_and_return_greetings, message)
            print greeting

    @staticmethod
    def _create_and_return_greetings(tx, message):
        result = tx.run("CREATE (a:Greeting) "
                        "SET a.message = $message "
                        "RETURN a.message +', from node ' + id(a)", message=message)
        return result.single()[0]

if __name__ == "__main__":
    neo4j_instance = HelloWorldExample('bolt://localhost:7687', 'neo4j', 'sergvcx')
    neo4j_instance.print_greetings('Hello World!')

