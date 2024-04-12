#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-12 at 09:59

@author: cook
"""

from typing import Any, Dict

from apero_checks import raw_tests
from apero_checks import red_tests
from apero_checks.core import general
from apero_checks.core import base
from apero_checks.core import io
from apero_checks.core import misc

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def find_override_test(params: Dict[str, Any]):
    # get test_name
    test_name = params['test_name']
    # test raw tests
    if test_name in raw_tests.override_list:
        return True, 'raw'
    # test red tests
    elif test_name in red_tests.override_list:
        return True, 'red'
    # otherwise do not allow override
    else:
        return False, None


def override_tests(params: Dict[str, Any],
                   test_type: str) -> Dict[str, Dict[str, Any]]:
    # get observation directories
    obsdirs = general.get_obs_dirs(params)
    # storage for test values
    test_values = dict()
    # get current data
    current_dataframe = general.get_current_dataframe(params, test_type)
    # index by obsdir
    obsdir_dataframe = current_dataframe.set_index('obsdir')
    # get test_name
    test_name = params['test_name']
    # -------------------------------------------------------------------------
    # loop around observation directories
    for obsdir in obsdirs:
        # print message on which observation directory we are processing
        msg = '*' * 50
        msg += f'\nProcessing observation directory {obsdir}'
        msg += '\n' + '*' * 50
        misc.log_msg(msg, level='info')
        # get a list of tests
        if test_type == 'raw':
            test_list = raw_tests.test_dict
        elif test_type == 'red':
            test_list = red_tests.test_dict
        else:
            emsg = 'RUN_TESTS error: test_type must be set to "raw" or "red"'
            raise base.AperoChecksError(emsg)
        # storage for this observation directory
        test_values[obsdir] = dict()
        # loop around all tests
        for it, test_name in enumerate(test_list):
            # deal with skipping tests
            if params['test_filter'] is not None:
                if test_name not in params['test_filter']:
                    if test_name in obsdir_dataframe:
                        value = obsdir_dataframe.loc[obsdir, test_name]
                    else:
                        value = ''
                    # push into test_values
                    test_values[obsdir][test_name] = value
        # ---------------------------------------------------------------------
        # Ask the user if they want to change the test value
        if test_name in test_values[obsdir]:
            # get the current value
            current_value = test_values[obsdir][test_name]
            # deal with non-boolean values
            if str(current_value) in ['', 'None', 'Null']:
                continue
            elif str(current_value).lower() in ['true', '1']:
                current_value = True
            else:
                current_value = False
            # construct the question to ask user
            question = (f'Current value for test {test_name} is '
                        f'{current_value}. Change to {not current_value}?'
                        f'[Y]es or [N]o\t')
            # ask user
            uinput = str(input(question))
            # deal with user input
            if 'y' in uinput.lower():
                new_value = not current_value
            else:
                new_value = current_value
            # push new value into test_values
            test_values[obsdir][test_name] = new_value
        # ---------------------------------------------------------------------
    # return the test values
    return test_values


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
