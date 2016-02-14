"""
This module contains all the routines that are respojnsible for pulling
the matrixes out of the neo4j graph and processing them.
"""
import hashlib
import itertools
import json
import os
import pickle
import string
from collections import defaultdict
from copy import copy
from random import shuffle, sample
from time import time
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

from bioflow.utils.gdfExportInterface import GdfExportInterface
from bioflow.utils.io_routines import write_to_csv, dump_object, undump_object
from bioflow.utils.log_behavior import get_logger
from bioflow.main_configs import Chromosome_source, Chromosome_file_filter
from bioflow.main_configs import Dumps, Outputs, interactome_rand_samp
from bioflow.internal_configs import edge_type_filters, Adjacency_Martix_Dict, \
    Conductance_Matrix_Dict
from bioflow.neo4j_db.GraphDeclarator import DatabaseGraph
from bioflow.neo4j_db.db_io_routines import expand_from_seed, \
    erase_custom_fields, node_extend_once, get_bulbs_id
from bioflow.neo4j_db.graph_content import bulbs_names_dict
from bioflow.algorithms_bank import conduction_routines as cr
from bioflow.neo4j_db.db_io_routines import stable_get_all

log = get_logger(__name__)

# Debug Log:
#   Main connex set is 24k nodes, with 4331 UP links and 1051 Hint links
#   With the full UP import, the main connex set is 25k Nodes, with 4336
#   UP Links and 1191 HiNT links
# With forbidding overloaded items: 24k771 Nodes, 4293 UP Links, 1186 HiNT
# links


# TODO: CRITICAL proper forbidden Ids Integration in the reaction set

class InteractomeInterface(object):
    """
    Interface between interactome in the database and the interactome graph laplacian

    :param main_connex_only: if True, loads only the elements of the main connex
    :param full_impact: if True, will compare the names of entities and link the
    entities that has the name with the "Possibily the same" relation
    """

    # List of all the reaction types present in the DatabaseGraph that will be
    # used as roots to build the interaction network (not all nodes are necessary
    # within the connex part of the graph)

    reactions_types_list = ['TemplateReaction', 'Degradation', 'BiochemicalReaction']
    reactions_types_list = [bulbs_names_dict[short_name][0] for short_name in reactions_types_list]

    def __init__(self, main_connex_only, full_impact):
        self.connexity_aware = main_connex_only
        self.full_impact = full_impact

        self.init_time = time()
        self.partial_time = time()

        self.adjacency_Matrix = np.zeros((4, 4))
        # This is just non-normalized laplacian matrix
        self.laplacian_matrix = np.zeros((4, 4))

        self.adj_eigenvects = np.zeros((4, 4))
        self.adj_eigenvals = np.zeros((4, 4))
        self.cond_eigenvects = np.zeros((4, 4))
        self.cond_eigenvals = np.zeros((4, 4))

        self.bulbs_id_2_matrix_index = {}
        self.matrix_index_2_bulbs_id = {}
        self.bulbs_id_2_display_name = {}
        self.bulbs_id_2_legacy_id = {}
        self.bulbs_id2_node_type = {}
        self.bulbs_id_2_localization = {}
        self.reached_uniprots_bulbs_id_list = []
        self.all_uniprots_bulbs_id_list = []
        self.Uniprot_attachments = {}  # currently maintained for legacy reasons
        self.uniprot_matrix_index_list = []

        # the properties below are used for iterative expansion of the dictionary
        self.ReactLinks = []
        self.InitSet = []
        self.GroupLinks = {}
        self.GroupSet = []
        self.sec_links = {}
        self.sec_set = []
        self.UP_Links = {}
        self.UPSet = []
        self.pre_full_set = []
        self.hint_links = {}
        self.full_set = []
        self.Super_Links = {}
        self.biogrid_links = {}
        self.ExpSet = []
        self.Highest_Set = []

        self.maps_dumps_location = Dumps.interactome_maps
        self.adjacency_dumps_location = Dumps.interactome_adjacency_matrix
        self.laplacian_dumps_location = Dumps.interactome_laplacian_matrix

        char_set = string.ascii_uppercase + string.digits
        self.random_tag = ''.join(sample(char_set * 6, 6))

        self.entry_point_uniprots_bulbs_ids = []
        self.UP2UP_voltages = {}
        self.uniprots_2_voltage_and_circulation = {}  # does not really seem to be used
        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}

        self.UP2Chrom = {}
        self.chromosomes_2_uniprot = defaultdict(list)

        self.incomplete_compute = False  # used in case of sparse sampling
        self.background = []
        self.connected_uniprots = []

        self.main_set = self.bulbs_id_2_matrix_index

    def pretty_time(self):
        """
        Times the execution

        :return: tuple containing the time since the creation of the Matrix_getter object and
         since the last cal of function formatted as string
        """
        it, pt = (round(time() - self.init_time),
                  round(time() - self.partial_time))
        pload = 'total: %s m %s s, \t partial: %s m %s s' % (
            int(it) / 60, it % 60, int(pt) / 60, pt % 60)
        self.partial_time = time()
        return pload

    def _time(self):
        """
        time since previous partial timer flag was reset.
        """
        pt = time() - self.partial_time
        return pt

    def dump_matrices(self):
        """
        dumps self.adjacency_Matrix and self.laplacian_matrix
        """
        dump_object(Dumps.interactome_adjacency_matrix, self.adjacency_Matrix)
        dump_object(Dumps.interactome_laplacian_matrix, self.laplacian_matrix)

    def undump_matrices(self):
        """
        undumps self.adjacency_Matrix and self.laplacian_matrix
        """
        self.adjacency_Matrix = undump_object(Dumps.interactome_adjacency_matrix)
        self.laplacian_matrix = undump_object(Dumps.interactome_laplacian_matrix)

    def dump_eigen(self):
        """
        dumps self.adj_eigenvals and self.laplacian_matrix and writes them to csv
        """
        write_to_csv(Dumps.eigen_VaMat, self.adj_eigenvals)
        write_to_csv(Dumps.eigen_ConMat, self.cond_eigenvals)
        dump_object(Dumps.val_eigen, (self.adj_eigenvals, self.adj_eigenvects))
        dump_object(
            Dumps.cond_eigen,
            (self.cond_eigenvals,
             self.cond_eigenvects))

    def undump_eigen(self):
        """
        undumps self.adj_eigenvals and self.laplacian_matrix
        """
        self.adj_eigenvals, self.adj_eigenvects = undump_object(
            Dumps.val_eigen)
        self.cond_eigenvals, self.cond_eigenvects = undump_object(
            Dumps.cond_eigen)

    def dump_maps(self):
        """
        dumps all the elements required for the mapping between the types and ids
         of database entries and matrix columns
        """
        log.debug("pre-dump e_p_u_b_i length: %s", len(self.entry_point_uniprots_bulbs_ids))
        log.debug("dumping into: %s", Dumps.interactome_maps)
        dump_object(
            Dumps.interactome_maps,
            (self.bulbs_id_2_matrix_index,
             self.matrix_index_2_bulbs_id,
             self.bulbs_id_2_display_name,
             self.bulbs_id_2_legacy_id,
             self.bulbs_id2_node_type,
             self.bulbs_id_2_localization,
             self.reached_uniprots_bulbs_id_list,
             self.all_uniprots_bulbs_id_list,
             self.Uniprot_attachments,
             self.UP2Chrom,
             self.chromosomes_2_uniprot,
             self.uniprot_matrix_index_list,
             self.entry_point_uniprots_bulbs_ids))  # TODO: delete here and below entry point u_b_i

    def undump_maps(self):
        """
        undumps all the elements required for the mapping between the types and ids of
        database entries and matrix columns
        """
        log.debug("undumping from %s", Dumps.interactome_maps)
        self.bulbs_id_2_matrix_index, self.matrix_index_2_bulbs_id, \
            self.bulbs_id_2_display_name, self.bulbs_id_2_legacy_id, self.bulbs_id2_node_type, \
            self.bulbs_id_2_localization, self.reached_uniprots_bulbs_id_list,\
            self.all_uniprots_bulbs_id_list, self.Uniprot_attachments, self.UP2Chrom, \
            self.chromosomes_2_uniprot, self.uniprot_matrix_index_list,\
            self.entry_point_uniprots_bulbs_ids = \
            undump_object(Dumps.interactome_maps)
        log.debug("post-undump e_p_u_b_i length: %s", len(self.entry_point_uniprots_bulbs_ids))

    def dump_memoized(self):
        md5 = hashlib.md5(
            json.dumps(
                sorted(
                    self.entry_point_uniprots_bulbs_ids),
                sort_keys=True)).hexdigest()

        payload = {
            'UP_hash': md5, 'sys_hash': self.md5_hash(),
            'size': len(self.entry_point_uniprots_bulbs_ids),
            'UPs': pickle.dumps(self.entry_point_uniprots_bulbs_ids),
            'currents': pickle.dumps((self.current_accumulator, self.node_current)),
            'voltages': pickle.dumps(self.uniprots_2_voltage_and_circulation)}
        dump_object(Dumps.Interactome_Analysis_memoized, payload)

    @staticmethod
    def undump_memoized():
        """ undumps memoized analysis """
        return undump_object(Dumps.Interactome_Analysis_memoized)

    # TODO: we can refactor that element as a separate object, responsible only for
    # laplacian_pulling
    # in it's current state, the complexity of this is too high (17)
    def full_load_ls(self):
        """
        Performs the loading of all the objects and the relations in the database
         according to the following algorithm:
            - get all the reaction nodes, i.e. reactions specified in the self.Reaction_List.
            - for each reaction, add the participants to the "seeds" list and create a set of
            groups linking them all together
            - for each element of hte "seeds" list, crawl the objects reachable from this seed
            by specific relations, add those new elements to the list and remember as reachable
            from the current element of the seed list.

            The relations are crawled in the following order:
                * "Group"
                * "Same"
                * "Contact_interaction" (repeated 5 times to increase as much as possible
                the main connex size, even for the knowledge groups that are badly
                represented in the interactome database)
                * "HiNT Contact_interaction"
                * "possibly_same"
            The exact relations in each group of relation abbreviations are specified in the
            configs.py file,

        This function fills in the self.Set and self.Links elements of the MatrixGetter object
        """

        def get_reaction_blocks():
            """
            Recovers the blocks if interaction that are due to a common set of reactions
            for the elements. They will be used as roots to build the complete interaction
            tree later on.

            :return:* Reagent clusters (list of lists, where each list represents a set of
                    reagent bubls_ids participating to the same reaction),

                    * seeds list,

                    * number of elements in the seeds list
            """
            reagent_clusters = []
            seeds = set()
            total_reaction_participants = 0
            for ReactionType in self.reactions_types_list:
                if not stable_get_all(ReactionType):
                    continue
                for Reaction in stable_get_all(ReactionType):
                    if Reaction is None:
                        continue
                    reaction_participants, reaction_particpants_no = node_extend_once(
                        edge_type_filters["Reaction"], self.connexity_aware, Reaction)
                    total_reaction_participants += reaction_particpants_no

                    if len(reaction_participants) > 1:
                        reagent_clusters.append(copy(reaction_participants))
                        seeds.update(reaction_participants)

            return reagent_clusters, seeds, total_reaction_participants

        def get_expansion(sub_seed, edge_type_filter):
            """
            Recovers all the nodes reached from the sub_seed according to a the relations listed
            in edge_type_filter

            :param sub_seed: set of NodeIDs that serve as a root for searching further relations

            :param edge_type_filter: type of relations according to which the graph will
             be explored

            :return clusters: * clusters (Dictionary of lists, where key is the
                              NodeID from the sub_seed)

                              * super_seed (List of NodeIDs attained from the SubSet by
                              the relation types from the edge_type_filter)

                              * total_expansion_participants (number of elements added to
                              the super_seed compared to the sub_seed)
            """
            clusters = {}
            super_seed = set()
            super_seed.update(sub_seed)
            total_expansion_participants = 0
            for element in sub_seed:
                # TODO: refactor below: edge_type_filter -> edge_type_filters[edge_type]
                local_list, count_increase = expand_from_seed(
                    element, edge_type_filter, self.connexity_aware)
                if len(local_list) > 0:
                    clusters[element] = copy(local_list)
                    super_seed.update(local_list)
                    total_expansion_participants += count_increase
            return clusters, super_seed, total_expansion_participants

        def print_characteristics(name, links, group, total_expansion_participants):
            """
            Prints a charectristics of a return object from the get_reaction_blocks or get_expansion

            :param name: name of the object characterised
            :param links: element 1 of the return tuple (links between objects)
            :param group: element 2 of the return tuple (super_seed)
            :param total_expansion_participants: element 3 of the return tuple
            (total_expansion_participants)
            """
            log.info('===========> %s <===========', name)
            log.info("links: %s, group: %s, nodes reached in expansion: %s",
                     len(links), len(group), total_expansion_participants)
            log.info(self.pretty_time())

        def n_expansion(starting_group, edge_type, char_name, expansions_n=1):
            """
            performs a single round of expansion/characterization

            :param starting_group:
            :param edge_type:
            :param char_name:
            :param expansions_n:
            :return:
            """
            tmp_links, tmp_group = ({}, starting_group.copy())

            for _i in range(0, expansions_n):
                round_name = '%s %s' % (char_name, _i+1)
                tmp_links, tmp_group, _count = get_expansion(
                    tmp_group, edge_type_filters[edge_type])
                if expansions_n > 1:
                    print_characteristics(round_name, tmp_links, tmp_group, _count)
                else:
                    print_characteristics(char_name, tmp_links, tmp_group, _count)

            return tmp_links, tmp_group

        ##################################################################################

        log.info('Entering retrieval of the connection system of physical entities')

        self.ReactLinks, self.InitSet, count = get_reaction_blocks()
        print_characteristics('Reactions', self.ReactLinks, self.InitSet, count)

        self.GroupLinks, self.GroupSet = n_expansion(self.InitSet, 'Group', 'Groups')

        self.sec_links, self.sec_set = n_expansion(self.GroupSet, 'Contact_interaction',
                                                   'Secondary Links', 6)

        self.UP_Links, self.UPSet = n_expansion(self.sec_set, 'Same', 'Uniprot Links')

        self.hint_links, self.pre_full_set = n_expansion(self.UPSet, 'HiNT_Contact_interaction',
                                                         'HiNT Links', 5)

        self.biogrid_links, self.full_set = n_expansion(self.pre_full_set,
                                                        'BioGRID_Contact_interaction',
                                                        'biogrid_path Links', 2)

        self.Super_Links, self.ExpSet = n_expansion(self.full_set, 'possibly_same',
                                                    'Looks_similar Links')

    # TODO: complexity too high (11); needs to be reduced.
    # Easy: this method actually contains two methods.
    #   1) One that pulls the uniprots reachable within the reactome routes
    #   2) Second one that pulls all the uniprots from
    #
    def map_rows_to_names(self):
        """
        Maps Node Database IDs, Legacy IDs, display names and types to matrix row/column indexes;
        """

        def request_location(_location_buffer_dict, location):
            """
            Just a Buffered lookup of location, since the number of cellular location
            is relatively small (~80), it makes sense to buffer the IOs on it.
            Normally should be moved out as a buffering decorator

            :param _location_buffer_dict: Buffered location
            :param location: location Node Legacy ID we are willing to verify
            :return: displayName of the requested location
            """
            location = str(location)
            if location in _location_buffer_dict.keys():
                return _location_buffer_dict[location]
            else:
                generator = DatabaseGraph.Location.index.lookup(ID=location)
                if generator is not None:
                    for elt in generator:
                        _location_buffer_dict[location] = str(elt.displayName)
                        return str(elt.displayName)

        #######################################################################

        counter = 0
        location_buffer_dict = {}

        self.bulbs_id_2_matrix_index = {}
        self.matrix_index_2_bulbs_id = {}

        log.info('nodes in Highest Set: %s', len(self.Highest_Set))
        for bulbs_node_id in self.Highest_Set:
            self.bulbs_id_2_matrix_index[bulbs_node_id] = counter
            self.matrix_index_2_bulbs_id[counter] = bulbs_node_id
            node = DatabaseGraph.vertices.get(bulbs_node_id)
            self.bulbs_id_2_display_name[bulbs_node_id] = node.displayName
            self.bulbs_id2_node_type[bulbs_node_id] = node.element_type
            self.bulbs_id_2_legacy_id[bulbs_node_id] = node.ID
            if node.element_type == "UNIPROT":
                # TODO: there is a problem: there is no UNIPROT nodes reached during the expansion.
                self.reached_uniprots_bulbs_id_list.append(bulbs_node_id)
                self.uniprot_matrix_index_list.append(counter)
            if node.localization is not None:
                self.bulbs_id_2_localization[bulbs_node_id] = request_location(
                    location_buffer_dict, node.localization)
            counter += 1

        self.all_uniprots_bulbs_id_list += self.reached_uniprots_bulbs_id_list
        self.reached_uniprots_bulbs_id_list = list(set(self.reached_uniprots_bulbs_id_list))
        log.info("reached uniprots: %s", len(self.reached_uniprots_bulbs_id_list))

        up_generator = stable_get_all(DatabaseGraph.UNIPORT)
        if up_generator:
            for up_node in up_generator:
                bulbs_node_id = get_bulbs_id(up_node)
                if bulbs_node_id not in self.reached_uniprots_bulbs_id_list:
                    self.all_uniprots_bulbs_id_list.append(bulbs_node_id)
                    self.bulbs_id_2_display_name[bulbs_node_id] = up_node.displayName
                    self.bulbs_id2_node_type[bulbs_node_id] = up_node.element_type
                    self.bulbs_id_2_legacy_id[bulbs_node_id] = up_node.ID

        self.all_uniprots_bulbs_id_list = list(set(self.all_uniprots_bulbs_id_list))
        log.info("analyzable uniprots: %s", len(self.all_uniprots_bulbs_id_list))

    def fast_row_insert(self, element, index_type):
        """
        Performs an correct insertion of an edge to the matrix.

        :param element: tuple of indexes designating elements we are willing to link
        :param index_type: type of the insert, so that the matrix coefficient can be
        looked up in the Adjacency_Martix_Dict or Conductance_Martix_Dict from the configs file
        """
        self.adjacency_Matrix[element[0], element[1]] =\
            min(self.adjacency_Matrix[element[0], element[1]] +
                Adjacency_Martix_Dict[index_type],
                1)

        self.adjacency_Matrix[element[1], element[0]] = \
            min(self.adjacency_Matrix[element[1], element[0]] +
                Adjacency_Martix_Dict[index_type],
                1)

        self.laplacian_matrix[element[0], element[1]] -= \
            Conductance_Matrix_Dict[index_type]
        self.laplacian_matrix[element[1], element[0]] -= \
            Conductance_Matrix_Dict[index_type]
        self.laplacian_matrix[element[1], element[1]] += \
            Conductance_Matrix_Dict[index_type]
        self.laplacian_matrix[element[0], element[0]] += \
            Conductance_Matrix_Dict[index_type]

    def create_val_matrix(self):
        """
        Creates anew self.adjacency_Matrix and self.laplacian_matrix
        """

        def insert_expansion_links(link_dict, link_type):
            for _key, values in link_dict.iteritems():
                for _val in values:
                    _link_indexes = (
                        self.bulbs_id_2_matrix_index[_key],
                        self.bulbs_id_2_matrix_index[_val])
                    self.fast_row_insert(_link_indexes, link_type)

        self.full_load_ls()

        self.Highest_Set = self.full_set

        if self.full_impact:
            self.Highest_Set = self.ExpSet

        self.map_rows_to_names()

        log.info("mapped bulbs_ids_to_ %s", self.pretty_time())

        log.info("building the ValMatrix %s ", self.pretty_time())

        load_len = len(self.Highest_Set)
        self.adjacency_Matrix = lil_matrix((load_len, load_len))
        self.laplacian_matrix = lil_matrix((load_len, load_len))

        for reaction_participant_set in self.ReactLinks:
            for reaction_participant in itertools.combinations(reaction_participant_set, 2):
                link_indexes = (
                    self.bulbs_id_2_matrix_index[reaction_participant[0]],
                    self.bulbs_id_2_matrix_index[reaction_participant[1]])
                self.fast_row_insert(link_indexes, "Reaction")

        insert_expansion_links(self.GroupLinks, "Group")
        insert_expansion_links(self.sec_links, "Contact_interaction")
        insert_expansion_links(self.UP_Links, "Same")
        insert_expansion_links(self.hint_links, "Contact_interaction")
        insert_expansion_links(self.biogrid_links, "weak_contact")
        insert_expansion_links(self.biogrid_links, "possibly_same")

    def get_eigen_spectrum(self, biggest_eigvals_to_get):
        """
        Recovers the eigenspectrum associated to the *n* biggest eigenvalues, where *n* is
        specified by biggest_eigvals_to_get. If the Adjacency and conductance matrix haven't
        been preloaded first, will raise an Exception

        :param biggest_eigvals_to_get: specifies how many biggest eigenvalues we are willing to get.
        :raise Exception: "Matrix must be pre-loaded first" if self.adjacency_Matrix and
        self.laplacian_matrix have not been computed anew or pre-loaded first
        """
        if self.laplacian_matrix.shape == (4, 4):
            log.critical("Matrix must be pre-loaded first")
            raise Exception("Matrix must be pre-loaded first")

        log.info("entering eigenvect computation; %s", self.pretty_time())

        self.adj_eigenvals, self.adj_eigenvects = eigsh(
            self.adjacency_Matrix, biggest_eigvals_to_get)
        self.cond_eigenvals, self.cond_eigenvects = eigsh(
            self.laplacian_matrix, biggest_eigvals_to_get)

        log.debug("Adjacency matrix eigenvalues:")
        log.debug(self.adj_eigenvals)
        log.debug('<======================>')
        log.debug("Laplacian matrix eigenvalues:")
        log.debug(self.cond_eigenvals)
        log.debug('<======================>')
        log.debug("all laplacian eigenvalues above%s", np.all(eigsh(self.laplacian_matrix)[0] > 0))
        log.info("Finished eigenvalues computation, starting the dump %s", self.pretty_time())

    def get_normalized_laplacian(self, fudge=1E-18):
        """
        Returns the normalized conductance matrix (Conductance matrix is actually the
         matrix laplacian of the associated Graph)

        :param fudge: Cancels out the action of eigenvalues that are too low.
        :return: normalized matrix
        """
        eigenvalues, eigenvectors = np.linalg.eigh(self.laplacian_matrix)
        diagonal = np.diag(1. / np.sqrt(eigenvalues + fudge))
        return np.dot(np.dot(eigenvectors, diagonal), eigenvectors.T)

    def write_connexity_infos(self):
        """
            Writes the infos about connexity of different components of a graph into the database
            for the future use. This execution is the main reason for the existance of the
            Adjacency Matrix.

            :warning: This process has to ber re-run each time the underlying database
            is changed in a way that migh affect the main connex graph
        """
        adjacency_graph_connected_components = connected_components(self.adjacency_Matrix,
                                                                    directed=False)

        counters = np.zeros((adjacency_graph_connected_components[0], 1))

        for elt in range(0, len(adjacency_graph_connected_components[1])):
            counters[adjacency_graph_connected_components[1][elt], 0] += 1

        giant_connex_component_index = np.argmax(counters)
        markings_to_do = len(adjacency_graph_connected_components[1])/10.

        for i, connex_component_index in enumerate(adjacency_graph_connected_components[1]):

            if i % int(markings_to_do) == 0:
                log.info('Marking graph main connex elements: %s done',
                         str("{0:.2f}".format(i / markings_to_do * 10.)))

            if connex_component_index == giant_connex_component_index:
                bulbs_node = DatabaseGraph.vertices.get(self.matrix_index_2_bulbs_id[i])
                bulbs_node.custom = 'Main_Connex'
                bulbs_node.main_connex = True
                bulbs_node.save()

        log.info("Marking of %s nodes for connexity was done in %s",
                 str(markings_to_do), str(self.pretty_time()))

    def full_rebuild(self):
        """
        Performs the initial loading routines that set up in place the sytem of dump files
        to allow fast loading in the future
        """
        connexity_aware = self.connexity_aware
        self.connexity_aware = False

        self.create_val_matrix()

        erase_custom_fields()
        self.write_connexity_infos()

        self.compute_uniprot_attachments()
        self.connexity_aware = connexity_aware

        self.create_val_matrix()
        self.get_eigen_spectrum(100)

        # self.map_UPs_to_chromosomes()  #TODO: revert once the chromosome
        # downloading from BioMart and assembly has been succesfully completed

        self.dump_maps()
        self.dump_matrices()
        self.dump_eigen()

    def fast_load(self):
        """
        Pre-loads mappings, matrices and eigenvalues
        """
        self.undump_maps()
        self.undump_matrices()
        self.undump_eigen()

        self.connected_uniprots = [
            NodeID for NodeID,
            idx in self.bulbs_id_2_matrix_index.iteritems() if idx < (
                self.laplacian_matrix.shape[0] - 1)]

    def get_descriptor_for_index(self, index):
        """
        Recovers a descriptor set for a given index in the current matrix mapping

        :param index: idenx for which return a descirptor
        :return: Type, displayName and if a localization is given, returns display name too.
        :rtype: tuple
        """
        if self.matrix_index_2_bulbs_id[index] in self.bulbs_id_2_localization.keys():
            return (self.bulbs_id2_node_type[self.matrix_index_2_bulbs_id[index]],
                    self.bulbs_id_2_display_name[self.matrix_index_2_bulbs_id[index]],
                    self.bulbs_id_2_localization[self.matrix_index_2_bulbs_id[index]])
        else:
            return (self.bulbs_id2_node_type[self.matrix_index_2_bulbs_id[index]],
                    self.bulbs_id_2_display_name[self.matrix_index_2_bulbs_id[index]])

    def compute_uniprot_attachments(self):
        """
        Computes the dictionary of attachments between the reached_uniprots_bulbs_id_list and
        Reactome proteins
        """
        log.info('attaching reactome proteins to uniprot nodes')
        uniprot_attachments_counter = 0
        reactome_attachments_counter = 0

        for uniprot_bulbs_id in self.reached_uniprots_bulbs_id_list:
            bulbs_node = DatabaseGraph.UNIPORT.get(uniprot_bulbs_id)
            reactome_bulbs_id_generator = bulbs_node.bothV("is_same")
            if reactome_bulbs_id_generator is not None:
                self.Uniprot_attachments[uniprot_bulbs_id] = []
                for uniprot_alias in reactome_bulbs_id_generator:
                    self.Uniprot_attachments[uniprot_bulbs_id].append(get_bulbs_id(uniprot_alias))
                uniprot_attachments_counter += 1
                reactome_attachments_counter += len(self.Uniprot_attachments[uniprot_bulbs_id])
                log.debug('attached %s Reactome proteins to the node %s',
                          len(self.Uniprot_attachments[uniprot_bulbs_id]), uniprot_bulbs_id)
            else:
                log.debug('No attachment for the node %s', uniprot_bulbs_id)

        log.info('Attached %s reactome protein nodes to %s / %s uniprot nodes',
                 reactome_attachments_counter,
                 uniprot_attachments_counter, len(self.reached_uniprots_bulbs_id_list))

    def hacky_corr(self):
        """
        Hacky method that should remain unused but by the devs.
        Generate the uniprot_matrix_index_list from the loaded Uniprot List.
        Prevents a 90-minute full database reload to compute something quickly.
        """
        self.undump_maps()
        self.uniprot_matrix_index_list = []
        for SP_Id in self.reached_uniprots_bulbs_id_list:
            bulbs_node = DatabaseGraph.UNIPORT.get(SP_Id)
            if bulbs_node.main_connex:
                self.uniprot_matrix_index_list.append(self.bulbs_id_2_matrix_index[SP_Id])
        log.info('number of indexed uniprots: %s', len(self.uniprot_matrix_index_list))
        self.dump_maps()

    def md5_hash(self):
        """
        Return the MD hash of self to ensure that all the defining properties have been correctly
        defined before dump/retrieval
        """
        sorted_initial_set = sorted(self.bulbs_id_2_matrix_index.keys())
        connected_ups = sorted(self.connected_uniprots)
        data = [
            self.connexity_aware,
            sorted_initial_set,
            self.full_impact,
            connected_ups]
        md5 = hashlib.md5(json.dumps(data, sort_keys=True)).hexdigest()

        return str(md5)

    def set_uniprot_source(self, uniprots):
        """
        Sets the reached_uniprots_bulbs_id_list on which the circulation computation routines will
        be performed by the otehr methods. Avoids passing as argument large lists of parameters.

        :param uniprots: List of node IDs of the uniprots on which we would like to
        perform current computations
        :raise Warning: if the uniprots were not present in the set of GOs for which
        we built the system or had no GO attached to them
        """
        if not set(uniprots) <= set(self.bulbs_id_2_matrix_index.keys()):
            log.warn('Following reached uniprots bulbs_ids were not retrieved upon the '
                     'circulation matrix construction: \n %s',
                     (set(uniprots) - set(self.bulbs_id_2_matrix_index.keys())))

        self.entry_point_uniprots_bulbs_ids = \
            [uniprot for uniprot in uniprots if uniprot in self.bulbs_id_2_matrix_index.keys()]

    # TODO: extract as the element performing a computation of the network
    # critical control parameters:
    #   - memoized => dismiss
    #   - sourced => dismiss
    #   - incremental => to be kept. So taht we can calculate the bacground even better
    #   - cancellation => to be kept
    #   - sparse_samples => to be kept
    #   - factor in the sampling run into this
    def build_extended_conduction_system(
            self,
            memoized=True,
            sourced=False,
            incremental=False,
            cancellation=True,
            sparse_samples=False):
        """
        Builds a conduction matrix that integrates uniprots, in order to allow an easier
        knowledge flow analysis

        :param bool memoized: if the tensions and individual relation matrices should be
            stored in the matrix and dumped at the end computation (required for submatrix
            re-computation)
        :param bool sourced: if true, all the relations will be looked up and not computed.
            Useful for the retrieval of sub-circulation group, but requires the
            uniprots_2_voltage_and_circulation to be pre-filled
        :param bool incremental: if True, all the circulation computation will be added to the
            existing ones. Useful for the computation of particularly big systems with
            intermediate dumps
        :param bool cancellation: divides the final current by number of bioflow-sink pairs
        :param int sparse_samples: if set to an integer the sampling will be sparse and not dense,
            i.e. instead of computation for each node pair, only an estimation will be made, equal to
            computing sparse_samples association with other randomly chosen nodes
        :return: adjusted conduction system
        """
        if not incremental or self.current_accumulator == np.zeros((2, 2)):
            self.current_accumulator = lil_matrix(self.laplacian_matrix.shape)
            self.UP2UP_voltages = {}
            self.node_current = defaultdict(float)
            if not sourced:
                self.uniprots_2_voltage_and_circulation = {}

        if sparse_samples:
            current_accumulator = cr.sample_group_edge_current(
                self.laplacian_matrix,
                [
                    self.bulbs_id_2_matrix_index[UP] for UP in self.entry_point_uniprots_bulbs_ids],
                re_samples=sparse_samples,
                cancellation=cancellation)
        else:
            current_accumulator, up_pair_2_voltage_current =\
                cr.group_edge_current_memoized(
                    self.laplacian_matrix,
                    [self.bulbs_id_2_matrix_index[UP]
                     for UP in self.entry_point_uniprots_bulbs_ids],
                    cancellation=cancellation,
                    # memoized=memoized,
                    memory_source=self.uniprots_2_voltage_and_circulation)
            # self.uniprots_2_voltage_and_circulation.update(up_pair_2_voltage_current)
            self.UP2UP_voltages.update(
                dict((pair, voltage)
                     for pair, (voltage, current) in up_pair_2_voltage_current.iteritems()))

        if incremental:
            self.current_accumulator = self.current_accumulator + current_accumulator
        else:
            self.current_accumulator = current_accumulator

        if memoized:
            self.dump_memoized()

        index_current = cr.get_current_through_nodes(self.current_accumulator)
        log.info('current accumulator shape %s', current_accumulator.shape)
        self.node_current.update(
            dict(
                (self.matrix_index_2_bulbs_id[idx],
                 val) for idx,
                val in enumerate(index_current)))

    def format_node_props(self, node_current, limit=0.01):
        """
        Formats the nodes for the analysis by in the knowledge_access_analysis module

        :param node_current: Current through the entity
        :param limit: hard limit to go_namespace_filter out the GO terms with too little current
        (compensates the minor currents in the gird)
        :return: {Entity:[node current, node degree]}
        """
        characterization_dict = {}
        limit_current = max(node_current.values()) * limit
        for NodeID, i in self.bulbs_id_2_matrix_index.iteritems():
            if node_current[NodeID] > limit_current:
                characterization_dict[NodeID] = [node_current[NodeID],
                                                 self.laplacian_matrix[i, i]]
        return characterization_dict

    def export_conduction_system(self, output_location=Outputs.Interactome_GDF_output):
        """
        Computes the conduction system of the GO terms and exports it to the GDF format and
         flushes it into a file that can be viewed with Gephi

         :param output_location:
        """

        if self.incomplete_compute:
            log.warning('Computation of the information circulation was not complete, %s',
                        'most likely due to the sampling')

        node_char_names = [
            'Current',
            'Type',
            'Legacy_ID',
            'Names',
            'Degree',
            'Source']

        node_char_types = [
            'DOUBLE',
            'VARCHAR',
            'VARCHAR',
            'VARCHAR',
            'DOUBLE',
            'DOUBLE']

        characterization_dict = {}

        log.info("Laplacian shape %s", self.laplacian_matrix.shape)
        log.info("Matrix size %s", max(self.matrix_index_2_bulbs_id.iterkeys()))

        for NodeID in self.node_current.iterkeys():
            matrix_index = self.bulbs_id_2_matrix_index[NodeID]
            characterization_dict[NodeID] = [
                str(self.node_current[NodeID]),
                self.bulbs_id2_node_type[NodeID],
                self.bulbs_id_2_legacy_id[NodeID],
                self.bulbs_id_2_display_name[NodeID].replace(',', '-'),
                str(self.laplacian_matrix[matrix_index, matrix_index]),
                str(float(int(NodeID in self.entry_point_uniprots_bulbs_ids)))]

        gdf_exporter = GdfExportInterface(
            target_fname=output_location,
            field_names=node_char_names,
            field_types=node_char_types,
            node_properties_dict=characterization_dict,
            min_current=0.01,
            index_2_label=self.matrix_index_2_bulbs_id,
            label_2_index=self.bulbs_id_2_matrix_index,
            current_matrix=self.current_accumulator)
        gdf_exporter.write()

    def export_subsystem(self, uniprot_system, uniprot_subsystem):
        """
        Exports the subsystem of reached_uniprots_bulbs_id_list and circulation between
         them based on a larger precalculated system.This is possible only of the memoization
         parameter was on during the execution of "build_extended_circulation_system()"
        function execution.

        :param uniprot_system: The set of uniprots for which the larger system was calculated
        :param uniprot_subsystem: the set of reached_uniprots_bulbs_id_list we are interested in
        :raise Exception: if the set of uniprots for which the larger system was calculated
         doesn't correspond to what is stored in the dumps
        """
        current_recombinator = self.undump_memoized()
        if not set(uniprot_system) == set(
                pickle.loads(current_recombinator['UPs'])):
            raise Exception('Wrong UP system re-analyzed')
        self.uniprots_2_voltage_and_circulation = pickle.loads(
            current_recombinator['voltages'])
        self.set_uniprot_source(uniprot_subsystem)
        self.build_extended_conduction_system(memoized=False, sourced=True)
        self.export_conduction_system()

    # TODO: remove memoization: it is not really used anywhere
    # parameters to remove:
    #   chromosome_specific (not now) we might want instead to implement it otherwise later
    #   memoized => we will never use ti
    def randomly_sample(
            self,
            samples_size,
            samples_each_size,
            sparse_rounds=False,
            chromosome_specific=False,
            memoized=False,
            no_add=False):
        """
        Randomly samples the set of reached_uniprots_bulbs_id_list used to create the model.
         This is the null model creation routine


        :param samples_size: list of numbers of uniprots we would like to create the model for
        :param samples_each_size: how many times we would like to sample each uniprot number
        :param sparse_rounds:  if we want to use sparse sampling (useful in case of
        large uniprot sets), we would use this option
        :param chromosome_specific: if we want the sampling to be chromosome-specific,
        set this parameter to the number of chromosome to sample from
        :param memoized: if set to True, the sampling would be remembered for export. Useful in
        case of the chromosome comparison
        :param no_add: if set to True, the result of sampling will not be added to the
        database of samples. Useful if re-running tests with similar parameters several times.
        :raise Exception: if the number of items in the samples size ann samples_each size
        are different
        """
        if not len(samples_size) == len(samples_each_size):
            raise Exception('Not the same list sizes!')

        if self.background:
            self.connected_uniprots = list(
                set(self.connected_uniprots).intersection(set(self.background)))

        if chromosome_specific:
            self.connected_uniprots = list(set(self.connected_uniprots).intersection(
                set(self.chromosomes_2_uniprot[str(chromosome_specific)])))

        for sample_size, iterations in zip(samples_size, samples_each_size):
            for i in range(0, iterations):
                shuffle(self.connected_uniprots)
                analytics_uniprot_list = self.connected_uniprots[:sample_size]
                self.set_uniprot_source(analytics_uniprot_list)
                self.build_extended_conduction_system(
                    memoized=memoized, sourced=False, sparse_samples=sparse_rounds)

                md5 = hashlib.md5(
                    json.dumps(
                        sorted(analytics_uniprot_list),
                        sort_keys=True)).hexdigest()

                if not no_add:
                    interactome_rand_samp.insert(
                        {
                            'UP_hash': md5,
                            'sys_hash': self.md5_hash(),
                            'size': sample_size,
                            'chrom': str(chromosome_specific),
                            'sparse_rounds': sparse_rounds,
                            'UPs': pickle.dumps(analytics_uniprot_list),
                            'currents': pickle.dumps(
                                (self.current_accumulator,
                                 self.node_current)),
                            'voltages': pickle.dumps(
                                self.UP2UP_voltages)})

                log.info('Random ID: %s \t Sample size: %s \t iteration: %s\t compops: %s \t '
                         'time: %s ', self.random_tag, sample_size, i,
                         "{0:.2f}".format(sample_size * (sample_size - 1) / 2 / self._time()),
                         self.pretty_time())

    # TODO: this method needs to be refactored in order to use the new chromosome to up id sources.
    # it needs to support the fact that chromosome mapping file are not necessarily present
    def map_uniprots_to_chromosomes(self):
        """
        Maps the Uniprot Ids to chromosomes of the organism and reversely.
        """
        resulting_dict = {}
        # TODO: add the possibility that Chromosome_source wasn't properly declared or is
        # a mapping to a void declaration if this is the case,
        for fle in os.listdir(Chromosome_source):
            if Chromosome_file_filter in fle:
                source_fle = Chromosome_source + '/' + fle
                chromosome_id = fle.split('.')[0].split(Chromosome_file_filter)[1]
                resulting_dict[chromosome_id] = open(source_fle).read()

        for i, (bulbs_id, legacy_id) in enumerate(self.bulbs_id_2_legacy_id.iteritems()):
            for chromosome_id, chromosome_text_block in resulting_dict.iteritems():
                if legacy_id in chromosome_text_block:
                    self.UP2Chrom[bulbs_id] = chromosome_id
                    self.chromosomes_2_uniprot[chromosome_id].append(bulbs_id)

        self.chromosomes_2_uniprot = dict(self.chromosomes_2_uniprot)


if __name__ == "__main__":
    interactome_interface_instance = InteractomeInterface(main_connex_only=True,
                                                          full_impact=False)
    # interactome_interface_instance.full_rebuild()
    interactome_interface_instance.fast_load()
    # print interactome_interface_instance.pretty_time()
    # print interactome_interface_instance.md5_hash()

    # interactome_interface_instance.set_Uniprot_source(test_set)
    # interactome_interface_instance.export_subsystem(test_set, test2)
    # interactome_interface_instance.build_extended_conduction_system()
    # interactome_interface_instance.export_conduction_system()
    # interactome_interface_instance.randomly_sample([100,250],[5,5], sparse_rounds=10)
