"""
Module containing an object that allows an easy export of the matrix-encoded information to GDF
without the need to export the whole relation matrix
"""
import numpy as np
from scipy.sparse import lil_matrix
from bioflow.configs.main_configs import Dumps
from bioflow.utils.general_utils.high_level_os_io import mkdir_recursive

# TODO: this class needs to be split into the GDF core that does all the work on a matrix rendering
# and filters that go_namespace_filter out the unwanted variables


class GdfExportInterface(object):
    """
    An interface allowing the export of the matrix relatioin object and node characteristics to a
    GDF format, compatible with visualization with appropriate tools.

    :param target_fname: name of the file to which the GDF file will be written to
    :param field_names: Names of different fields for the node description
    :param field_types: Types of different nodes for the node description
    :param node_properties_dict: dictionary mapping the node labels to outputs
    :param min_current: minimal current below which we are not rendering the links anymore
    :param index_2_label: Mapping from the indexes curent matrix lines/columns to the node labels
    :param current_matrix: matrix of currents from which we wish to rendred the GDF
    """

    Authorised_names = ['VARCHAR', 'DOUBLE', 'BOOLEAN']

    def __init__(
            self,
            target_fname,
            field_names,
            field_types,
            node_properties_dict,
            min_current,
            index_2_label,
            label_2_index,
            current_matrix,
            directed=False):
        mkdir_recursive(target_fname)
        self.target_file = open(target_fname, 'wt')
        self.field_types = field_types
        self.field_names = field_names
        self.node_properties = node_properties_dict
        self.Idx2Label = index_2_label
        self.Label2Idx = label_2_index
        self.current_Matrix = lil_matrix(current_matrix)
        # matrix where M[i,j] = current intesitu from i to j. Triangular superior, if current is
        #  from j to i, current is negative

        # current retrieval for the output should be done by getting all the non-zero terms of
        # the current matrix and then filtering out terms/lines that have too little absolute
        # current

        # rebuilding a new current Matrix and creating a dict to map the relations from the
        # previous matrix into a new one.
        self.mincurrent = min_current * \
            self.current_Matrix[self.current_Matrix.nonzero()].toarray().max()
        # minimal current for which we will be performing filtering out of the conductances and
        # nodes through which the traffic is below that limit
        self.directed = directed
        self.verify()

    def verify(self):
        """
        :raises Exception: "GDF Node declaration is wrong!" - if the length of field names and
        field type differ
        :raises Exception: "Wrong Types were declared, ...." -  if the declared types are not in
        the Authorised names
        """
        if len(self.field_types) != len(self.field_types):
            raise Exception('GDF Node declaration is wrong')
        if not set(self.Authorised_names) >= set(self.field_types):
            raise Exception(
                'Wrong types were declared. please refer to the Doc')

    def write_nodedefs(self):
        """
        Takes in the dictionary that maps and returns the nodedefs line

        """
        accumulator = []
        for node_name, node_type in zip(self.field_names, self.field_types):
            accumulator.append(node_name + ' ' + node_type)
        retstring = ', '.join(accumulator)
        retstring = 'nodedef> name VARCHAR, ' + retstring + '\n'
        self.target_file.write(retstring)

    def write_nodes(self):
        """
        Write the nodes with associated informations

        """
        for nodename, nodeprops in self.node_properties.items():
            if self.mincurrent and float(nodeprops[0]) > self.mincurrent:
                self.target_file.write(
                    str(nodename) + ', ' + ', '.join(nodeprops) + '\n')
            else:
                self.target_file.write(
                    str(nodename) + ', ' + ', '.join(nodeprops) + '\n')

    def write_edgedefs(self):
        """
        Write defintion for the edges. Right now, the information passing through the edges is
        restricted to the current

        """
        retstring = 'edgedef> node1 VARCHAR, node1 VARCHAR, weight DOUBLE, directed BOOLEAN\n'
        self.target_file.write(retstring)

    def write_edges(self):
        """
        Writes information about edges connections. This information are pulled from the
        conductance matrix.

        """
        nz = self.current_Matrix.nonzero()
        for i, j in zip(nz[0], nz[1]):
            if abs(self.current_Matrix[i, j]) > self.mincurrent:
                if self.directed:
                    write_line = ', '.join([str(self.Idx2Label[i]), str(self.Idx2Label[j]),
                                            str(self.current_Matrix[i, j]), 'true']) + '\n'
                else:
                    write_line = ', '.join([str(self.Idx2Label[i]), str(self.Idx2Label[j]),
                                            str(self.current_Matrix[i, j]), 'false']) + '\n'
                self.target_file.write(write_line)

    def write(self):
        """
        Performs all the writing routines and output file closing all at once

        """
        self.write_nodedefs()
        self.write_nodes()
        self.write_edgedefs()
        self.write_edges()
        self.target_file.close()


if __name__ == "__main__":
    pass
