from bioflow.annotation_network.knowledge_access_analysis import get_go_interface_instance
from bioflow.utils.linalg_routines import analyze_eigenvects


def linearly_independent_go_groups(size):
    """
    Performs the analysis of linearly independent GO groups

    :param size: number of most prominent eigenvectors contributing to an eigenvalue
    we would want to see
    """
    go_interface_instace = get_go_interface_instance()
    go_interface_instace._undump_independent_linear_sets()
    char_indexes = dict(
        (key,
         (len(
             go_interface_instace._limiter_go_2_up_reachable_nodes[value]),
             go_interface_instace.neo4j_id_2_legacy_id[value],
             go_interface_instace.neo4j_id_2_display_name[value])) for key,
                                                                       value in go_interface_instace.mat_idx_2_note_id.iteritems())
    print(go_interface_instace.pretty_time())
    analyze_eigenvects(go_interface_instace.indep_lapl, size, char_indexes)
