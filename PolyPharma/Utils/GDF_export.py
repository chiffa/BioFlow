__author__ = 'ank'


Authorised_names = ['VARCHAR', 'DOUBLE', 'BOOLEAN', 'DOUBLE']

class GDF_export_Interface(object):
    """

    :param target_fname: name of the file to which the GDF file will be written to
    :param field_names: Names of different fields for the node description
    :param field_types: Types of different nodes for the node description
    :param node_properties_dict:
    :param directed_adjacency_matrix:
    :param mincurrent:
    :param Label2idx:
    :param current_Matrix:
    """

    def __init__(self, target_fname, field_names, field_types, node_properties_dict, directed_adjacency_matrix, mincurrent, Label2idx, current_Matrix = None,):
        self.target_file = open(target_fname, 'w')
        self.field_types = field_types
        self.field_names = field_names
        self.node_properties = node_properties_dict
        self.directed_adjacency_matrix = node_properties_dict
        self.Label2idx = Label2idx
        self.mincurrent = mincurrent # minimal current for which we will be performing filtering out of the conductances and nodes through whichthe trafic is below that limit
        self.current_Matrix = current_Matrix # matrix where M[i,j] = current intesitu from i to j. Triangular superior, if current is from j to i, current is negative
        # current retrieval for the output should be done by getting all the non-zero terms of the current matrix and then filtering out terms/lines that have too little absolute current
        # rebuilding a new current Matrix and creating a dict to map the relations from the previous matrix into a new one.
        self.verify()


    def verify(self):
        """
        :raises Exception: "GDF Node declaration is wrong!" - of the length of field names and field type differ
        :raises Exception: "Wrong Types were declared, ...." - if the declared types are not in the Authorised names
        """
        if len(self.field_types) != len(self.field_types):
            raise Exception('GDF Node declaration is wrong')
        if not set(Authorised_names) >= set(self.field_types):
            raise  Exception('Wrong types were declared. please refer to the module for authorized types list')


    def write_nodedefs(self):
        """
        Takes in the dictionary that maps and returns the nodedefs line

        :returns: string line starting the nodedef line

        ..codeblock:: python

        """
        accumulator = []
        for name, type in zip(self.field_names, self.field_types):
            accumulator.append(name+' '+type)
        retstring = ', '.join(accumulator)
        retstring = 'nodedef> name VARCHAR, current DOUBLE, '+retstring+'\n'
        self.target_file.write(retstring)
        return retstring


    def write_nodes(self):
        pass


    def write_edgedefs(self):
        retstring = 'nodedef> node1 VARCHAR, node1 VARCHAR, current DOUBLE \n'
        self.target_file.write(retstring)
        return retstring


    def write_edges(self):
        pass


if __name__ == "__main__":
    GDFW = GDF_export_Interface('GDF_exporter_test.gdf', ['test'],['VARCHAR'], None, None, None, None)
    print GDFW.write_nodedefs()
