#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict

from apero_checks.core import io
from apero_checks.core import misc
from apero_checks.core import red_functions

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Test whether observation directory exists in the APERO raw data directory
    if it doesn't it means symlinks have not be created (i.e. trigger has not
    been run)

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # update apero profile
    apero_params = red_functions.update_apero_profile(params)
    # define parameters we use here
    raw_directory = apero_params['DRS_DATA_RAW']
    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------
    # get all observation directories
    obsdirs = io.get_obs_dirs(params)
    # -------------------------------------------------------------------------
    # test if observation directory exists in our list
    if obsdir not in obsdirs:
        if log:
            msg = ('OBSDIR TEST: Observation directory {0} does not exist in '
                   '{1} - TEST FAILED')
            margs = [obsdir, raw_directory]
            misc.log_msg(msg.format(*margs), level='warning')
        return False
    # -------------------------------------------------------------------------
    if log:
        msg = ('OBSDIR TEST: Observation directory {0} exists in {1} '
               '- TEST PASSED')
        margs = [obsdir, raw_directory]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------
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
    _obsdir = '2021-03-15'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
