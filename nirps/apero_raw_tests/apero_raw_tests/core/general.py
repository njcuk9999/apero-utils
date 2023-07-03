#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 17:03

@author: cook
"""
from typing import Any, Dict

from apero_raw_tests.tests import tests
from apero_raw_tests.core import core

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def run_test(params: Dict[str, Any], test_name: str, it: int, num_tests: int,
             log: bool = True):
    # try to run test
    try:
        # print which test we are running
        msg = 'Running test {0} [{1}/{2}]'
        margs = [test_name, it + 1, num_tests]
        core.log_msg(msg.format(*margs), level='')
        output = tests[test_name].test(params, log=log)
    except Exception as e:
        msg = '\tError in test {0} \n\t {1}'
        margs = [test_name, e]
        core.log_msg(msg.format(*margs), level='warning')
        output = False
    return output


def run_tests(params: Dict[str, Any]):
    """
    Run all tests in silent mode and return a dictionary of test values

    :param params:
    :return:
    """

    # storage for test values
    test_values = dict()

    # loop around all tests
    for it, test in enumerate(tests):
        # get the name of the test
        test_name = test.NAME
        # run the test
        output = run_test(params, test_name, it, len(tests), False)
        # add to test values
        test_values[test_name] = output
    # return the test values
    return test_values


def run_single_test(params: Dict[str, Any], test_name: str):
    """
    Run a single test in log mode

    :param params: Parameter dictionary containing constants
    :param test_name:  the name of the test to run
    :return:
    """
    if test_name in tests.tests(params):
        _ = run_test(params, test_name, 0, 1, True)


# -----------------------------------------------------------------------------
# Deal with uploading tests
def upload_tests(params: Dict[str, Any], results: Dict[str, bool]):
    pass


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
