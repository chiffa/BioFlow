"""
Contains a set of wrappers that are useful for debugging and operation profiling
"""
import matplotlib.pyplot as plt
from time import time
from BioFlow.utils.LogManager import logger


def render_2d_matrix(matrix, name):
    plt.title(name)
    plt.imshow(matrix, interpolation='nearest')
    plt.colorbar()
    plt.show()


def debug_wrapper(function_of_interest):

    def check_matrix(*args, **kwargs):
        result = function_of_interest(*args, **kwargs)
        if type(result) is not tuple:
            render_2d_matrix(result, function_of_interest.__name__)
        else:  # TODO: perform a vector rendering of function results
            render_2d_matrix(result[0], function_of_interest.__name__)
        check_matrix.__name__ = function_of_interest.__name__
        check_matrix.__doc__ = function_of_interest.__doc__
        return result

    return check_matrix


def time_it_wrapper(function_of_interest):

    def time_execution(*args, **kwargs):
        start = time()
        result = function_of_interest(*args, **kwargs)
        logger.debug('%s run in %s' % (function_of_interest.__name__, time() - start))
        time_execution.__name__ = function_of_interest.__name__
        time_execution.__doc__ = function_of_interest.__doc__
        return result

    return time_execution
