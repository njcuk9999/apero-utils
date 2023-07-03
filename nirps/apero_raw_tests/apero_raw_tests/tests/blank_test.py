#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict

# =============================================================================
# Define variables
# =============================================================================
# name of the column in database, must be unique
NAME = 'BLANK'

# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], log=False) -> bool:
    """
    Blank test - this tests whether test was run (should always return True)
    All other tests should return True or False, and only print messages if
    log is True.

    Paseed = True

    :param night:
    :param log:
    :return:
    """
    # define parameters from params (good practice)
    # must be in core.parameters.py
    night = params['obsdir']
    # -------------------------------------------------------------------------
    # all print out messages must be wrapped in if log
    if log:
        print('BLANK TEST: This is a blank test - it should always return True')
        print('NIGHT: {0}'.format(night))
    return True


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _params['obsdir'] = '2021-03-15'
    # run the test
    test(_params, log=True)

# =============================================================================
# End of code
# =============================================================================
