#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict

import os
import glob
from astropy.io import fits
from tqdm import tqdm

from apero_checks.core import apero_functions
from apero_checks.core import misc


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
    Test for pixel shifts in pp files (tmp directory) -
    Checks the DETOFFDX and DETOFFDY header keys to look for any non-zero values

    All 0 = True
    Any non-zero value = False

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """

    # update apero-profile
    apero_params = apero_functions.update_apero_profile(params)
    # tmp directory
    tmp_dir = apero_params['DRS_DATA_WORKING']
    # directory to check
    obsdir_path = os.path.join(tmp_dir, obsdir)

    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    # check if directory exists
    if not os.path.exists(obsdir_path):
        if log:
            print('tmp directory {} does not exist'.format(obsdir))
        return False

    # list of all the files in observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))

    # check if there are files in the directory
    if len(files) == 0:
        if log:
            print('No files in directory {}'.format(obsdir))
        return False

    passed = True
    failed_msg = ''

    # check pixel shift header keys for all files in obsdir
    for filename in tqdm(files, leave=False):

        hdr = fits.getheader(filename)
        filename = os.path.basename(filename)
        dx = float(hdr['DETOFFDX'])
        dy = float(hdr['DETOFFDY'])
        if dx != 0 or dy != 0:
            # there is a shift
            passed = False
            failed_msg += ('Shift dx = {}, dy = {} detected in file '
                           '{} \n').format(dx, dy, filename)

    if log:
        if len(failed_msg) > 0:
            print(failed_msg)
        else:
            print('No shifts detected')

    return passed


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
