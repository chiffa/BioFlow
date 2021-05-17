"""
Module containing functions that perform weighting policies
"""
from bioflow.utils.log_behavior import get_logger
import bioflow.configs.main_configs as confs

log = get_logger(__name__)


def flat_policy(start_node, end_node, edge) -> float:
    """
    A flat weight policy (any link is set to 1). Basically a connexity indicator

    :param start_node: start node object
    :param end_node: end node object
    :param edge: edge object
    :return: 1
    """
    return 1


def source_x_type_policy(start_node, end_node, edge, source_weights, type_weights) -> float:
    """
    A policy that uses the type and the source of the edge itself (and only the edge)

    :param start_node: start node object
    :param end_node: end node object
    :param edge: edge object
    :return: 1
    """
    return type_weights[edge.type] * source_weights[edge['source']]


def default_lapl_source_x_type_policy(start_node, end_node, edge) -> float:
    """
    A policy that uses the type and the source of the edge itself (and only the edge) and uses
    the laplacian dictionary weight maps hardcoded in the main_configs

    :param start_node: start node object
    :param end_node: end node object
    :param edge: edge object
    :return: weight
    """
    return source_x_type_policy(start_node, end_node, edge,
                                confs.laplacian_default_source_edge_weighting,
                                confs.laplacian_default_type_edge_weighting)


def default_adj_source_x_type_policy(start_node, end_node, edge) -> float:
    """
    A policy that uses the type and the source of the edge itself (and only the edge) and uses
    the laplacian dictionary weight maps hardcoded in the main_configs

    :param start_node: start node object
    :param end_node: end node object
    :param edge: edge object
    :return: weight
    """
    return source_x_type_policy(start_node, end_node, edge,
                                confs.adjacecency_default_source_edge_weighting,
                                confs.adjacency_default_type_edge_weighting)


active_default_lapl_weighting_policy = default_lapl_source_x_type_policy
active_default_adj_weighting_policy = default_adj_source_x_type_policy