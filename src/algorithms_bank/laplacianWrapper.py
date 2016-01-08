"""
Module containing laplacian routines wrapper that is used in combination with the conduction
routines to perform the analysis
"""
import pickle
import json
import hashlib
import numpy as np
from time import time
from itertools import izip

from src.utils.io_Routines import write_to_csv, dump_object, undump_object
from src.utils.log_behavior import get_logger

log = get_logger(__name__)


class LaplacianWrapper(object):

    def __init__(self):
        self.analytic_uniprots = []
        self.current_accumulator = np.zeros((2, 2))
        self.node_current = {}
        self.uniprots_2_voltage_and_circulation = {}

        self.memoization_dump_location = ''  # Dumps.Interactome_Analysis_memoized
        self.main_set = {}                   # sorted set that decides md_5 behavior
        self.md5_hash_variables = []         # properties of self that need to be in md5 hash

    def dump_memoized(self):
        """
        Dumps a memoization-compatible file
        :return:
        """
        md5 = hashlib.md5(
            json.dumps(
                sorted(
                    self.analytic_uniprots),
                sort_keys=True)).hexdigest()
        payload = {
            'UP_hash': md5, 'sys_hash': self.md5_hash(), 'size': len(
                self.analytic_uniprots), 'UPs': pickle.dumps(
                self.analytic_uniprots), 'currents': pickle.dumps(
                (self.current_accumulator, self.node_current)), 'voltages': pickle.dumps(
                    self.uniprots_2_voltage_and_circulation)}
        dump_object(self.memoization_dump_location, payload)

    def undump_memoized(self):
        """
        :return: undumped memoized analysis
        :rtype: dict
        """
        return undump_object(self.memoization_dump_location)

    def load_laplacian(self):
        """
        Loads the laplacian that was constructed
        """
        pass

    def dump_to_string(self):
        """
        Dumps the contents to the string from which the object can be rebuild
        """
        pass

    def undump_from_string(self):
        """
        Recovers the contents of the string and rebuilds the object from it.
        """
        pass

    def calculate_current_for_sample(self):
        """ calculates the current circulation for the sample """
        pass

    def sample_current(self):
        """
        Performs the current sampling for comparison.

        Conserves random node-pairs resistances and average inter-node tension in case of
        """
        pass

    def compare_circulation_to_reference(self):
        """
        Compares actual sample to a random sample from the background

        Compares the obtained current to the 95% estimation of chance that flow would be observed
        on a random sample of equal size analyzed with the same approach
        """
        pass

    def output_as_gdf(self):
        """ Returns the laplacian graph and it's circulation state """
        pass

    def retrieve_groups(self):
        """ Clusters together the nodes that are more critical than the """
        pass

