"""
Contains a set of wrappers that are useful for debugging and operation profiling
"""
import matplotlib.pyplot as plt
from time import time
from src.utils.log_behavior import get_logger

log = get_logger(__name__)


def render_2d_matrix(matrix, name):
    """
    Subroutine requried by the rendering wrapper.

    :param matrix:
    :param name:
    :return:
    """
    plt.title(name)
    plt.imshow(matrix, interpolation='nearest')
    plt.colorbar()
    plt.show()


def debug_wrapper(function_of_interest):
    """
    Convenient wrapper inspecting the results of a function returning a single matrix

    :param function_of_interest:
    :return: wrapped functions with copied name and documentation
    """

    def check_matrix(*args, **kwargs):
        result = function_of_interest(*args, **kwargs)
        if not isinstance(result, tuple):
            render_2d_matrix(result, function_of_interest.__name__)
        else:
            render_2d_matrix(result[0], function_of_interest.__name__)
        check_matrix.__name__ = function_of_interest.__name__
        check_matrix.__doc__ = function_of_interest.__doc__
        return result

    return check_matrix


def time_it_wrapper(function_of_interest):
    """
    Convenient wrapper for timing the execution time of functions
    :param function_of_interest:
    :return: wrapped functions with copied name and documentation
    """

    def time_execution(*args, **kwargs):
        start = time()
        result = function_of_interest(*args, **kwargs)
        log.debug('%s run in %s',
                  function_of_interest.__name__, time() - start)
        time_execution.__name__ = function_of_interest.__name__
        time_execution.__doc__ = function_of_interest.__doc__
        return result

    return time_execution


def my_timer(message='', previous_time=[]):
    """
    A small timer to be used in code to measure the execution duration of elements of code

    :param message:
    :param previous_time:
    :return:
    """
    if not previous_time:
        print 'set timer'
        previous_time.append(time())
    else:
        print '%s timer reset. Time since previous %s' % (message, time() - previous_time[0])
        previous_time[0] = time()
