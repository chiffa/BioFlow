from bioflow.configs import main_configs
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)


def check_primary_and_secondary_are_disjoined(primary_set, secondary_set):

    if secondary_set is None:
        return True

    if len(secondary_set) == 0:
        return True

    if set(primary_set).isdisjoint(set(secondary_set)):
        return True

    else:
        log.critical('Assumption on disjointness of the primary and secondary sets has failed.'
                     'Intersection database ids: %s ' % set(primary_set).intersection(set(secondary_set)))

        raise Exception('assumption failed')


def check_analytic_sets_assumptions(primary_set, secondary_set):

    # checks that fail will raise exceptions and interrupt the flow by itself
    check_primary_and_secondary_are_disjoined(primary_set, secondary_set)

    return True
