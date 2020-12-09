"""
This file contains scripts that are useful to analyze and profile the code behavior
"""

###################################################################################################
# line-by line performance profiler
import line_profiler
import atexit

profile = line_profiler.LineProfiler()
atexit.register(profile.print_stats)


@profile  # this is profiling wrapper
def function_to_test():
    # do something that is potentially a performance bottleneck
    pass

# after the execution the line-per line breakdown of the time each line takes will be printed to
# stdout
###################################################################################################


###################################################################################################
# external memory inspection
import os
import psutil

call_inc = 0
last_mem_log = 0

def log_mem(flag=''):
    if psutil_main_loop_memory_tracing:
        global last_mem_log
        global call_inc
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024 / 1024

        # print("DEBUG: memory load @ %s: %s" % (flag, mem))

        if call_inc == 0:
            print("DEBUG: memory load @ %s: %s" % (flag, mem))
        if mem != last_mem_log:
            print("DEBUG: memory load increased @ %s by %.2f on inc %d" % (flag, mem - last_mem_log,
                                                                           call_inc))
            last_mem_log = mem
    else:
        pass


def other_function_to_test():
    while True:  # loop that might leak memory
        log_mem('pre-location')
        # do something that is potentially a performance bottleneck
        log_mem('post-location')
        global call_inc
        call_inc += 1
    pass

# you will need to recover the output from the print and pass it through an external script to
# analyze where the memory is leaking.
# Beware of the location that seems to leak memory - it will be the location where the persistent
# object is created, not where it is retained. Because of that the line you might zero-in might
# not be the one causing problem.
###################################################################################################

###################################################################################################
# this will force all the warnings raised by python code to generate a traceback, so that you
# know where they come from and be raised every single time they arise rather than only the first
# time (default behavior)
import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


warnings.showwarning = warn_with_traceback
warnings.simplefilter("always")
###################################################################################################