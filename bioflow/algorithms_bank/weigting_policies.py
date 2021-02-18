from bioflow.utils.log_behavior import get_logger
import bioflow.configs.main_configs as confs

log = get_logger(__name__)


def flat_policy(start_node, end_node, edge) -> float:
    return 1


def source_x_type_policy(start_node, end_node, edge, source_weights, type_weights) -> float:
    return type_weights[edge.type] * source_weights[edge['source']]


def default_lapl_source_x_type_policy(start_node, end_node, edge) -> float:
    return source_x_type_policy(start_node, end_node, edge,
                                confs.laplacian_default_source_edge_weighting,
                                confs.laplacian_default_type_edge_weighting)


def default_adj_source_x_type_policy(start_node, end_node, edge) -> float:
    return source_x_type_policy(start_node, end_node, edge,
                                confs.adjacecency_default_source_edge_weighting,
                                confs.adjacency_default_type_edge_weighting)