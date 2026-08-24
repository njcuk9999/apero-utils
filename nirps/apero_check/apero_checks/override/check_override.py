#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-12 at 09:59

@author: cook
"""

from typing import Any, Dict, Tuple

from apero_checks import raw_tests
from apero_checks import red_tests
from apero_checks.core import general
from apero_checks.core import base
from apero_checks.core import misc

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
DictReturn = Dict[str, Dict[str, Any]]

def override_tests(params: Dict[str, Any],
                   test_type: str) -> Tuple[DictReturn, DictReturn]:
    # get observation directories
    obsdirs = general.get_obs_dirs(params)
    # storage for test values
    test_values = dict()
    # store overrides
    overrides = dict()
    # get current data
    current_dataframe = general.get_current_dataframe(params, test_type)
    # index by obsdir
    obsdir_dataframe = current_dataframe.set_index('obsdir')
    # get test_name
    test_name = params['test_name']
    # -------------------------------------------------------------------------
    # get any custom override values from the command line (None = ask user)
    override_current_value = params.get('current_value', None)
    override_who = params.get('who', None)
    override_comment = params.get('comment', None)
    # convert current_value from a string (True/False/None) to a bool or None
    if str(override_current_value) in ['None', 'Null', '']:
        override_current_value = None
    elif str(override_current_value).lower() in ['true', '1']:
        override_current_value = True
    else:
        override_current_value = False
    # if any of these were set, confirm with the user that they are correct
    if any(v is not None for v in [override_current_value, override_who,
                                   override_comment]):
        # build confirmation message
        cmsg = 'The following override values were provided on the command line:'
        cmsg += f'\n\t--current_value = {override_current_value}'
        cmsg += f'\n\t--who           = {override_who}'
        cmsg += f'\n\t--comment       = {override_comment}'
        cmsg += '\nIs this correct?\n [Y]es or [N]o\t'
        # ask user
        uinput = str(input(cmsg))
        # if not confirmed, reset to None so the questions are asked as normal
        if 'y' not in uinput.lower():
            override_current_value = None
            override_who = None
            override_comment = None
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
        # value
        overrides[obsdir] = dict()
        # loop around all tests
        for it, _test_name in enumerate(test_list):
            # get the test function
            if _test_name in obsdir_dataframe.columns:
                value = obsdir_dataframe.loc[obsdir, _test_name]
            else:
                value = ''
            # push into test_values
            test_values[obsdir][_test_name] = value
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
                        f'\n [Y]es or [N]o\t')
            # decide whether to apply the override
            #   if --current_value was given and it matches the actual current
            #   value, apply the override without asking, otherwise ask as normal
            if override_current_value is not None:
                do_override = override_current_value == current_value
                if not do_override:
                    do_override = 'y' in str(input(question)).lower()
            else:
                # ask user
                do_override = 'y' in str(input(question)).lower()
            # deal with user input
            if do_override:
                # get the name (from --who or ask the user)
                if override_who is not None:
                    name = override_who
                else:
                    name = ''
                    while len(name) == 0:
                        # ask for name
                        question = 'Please enter your name:\t'
                        name = str(input(question))
                # get the reason (from --comment or ask the user)
                if override_comment is not None:
                    reason = override_comment
                else:
                    reason = ''
                    while len(reason) == 0:
                        # ask for reason
                        question = 'Please enter the reason for the override:\t'
                        reason = str(input(question))
                # set the new value opposite to the old value
                new_value = not current_value
                # add to overrides
                overrides[obsdir][test_name] = (new_value, name, reason)
            else:
                new_value = current_value
            # push new value into test_values
            test_values[obsdir][test_name] = new_value
    # ---------------------------------------------------------------------
    # return the test values
    return test_values, overrides


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
