#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict, Tuple

import os
import glob
from astropy.io import fits
from tqdm import tqdm

from apero_checks.core import apero_functions
from apero_checks.core import misc
from apero_checks.core import io


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)
SNR_KEY = 'EXTSN064'
SNR_LIMIT = 10
SCI_DPRTYPES = ['OBJ_DARK', 'OBJ_FP', 'OBJ_SKY', 'TELLU_SKY', 'FLUXSTD_SKY']

# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> Tuple[bool, str]:
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
    red_dir = apero_params['DRS_DATA_REDUC']
    # directory to check
    obsdir_path = os.path.join(red_dir, obsdir)

    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    # check if directory exists
    if not os.path.exists(obsdir_path):
        out_msg = ('tmp directory {} does not exist'.format(obsdir))
        if log:
            print(out_msg)
        return False, out_msg

    # list of all the e2dsff A files in observation directory
    files = glob.glob(os.path.join(obsdir_path, '*pp_e2dsff_A.fits'))

    # check if there are files in the directory
    if len(files) == 0:
        out_msg = ('No files in directory {}'.format(obsdir))
        if log:
            print(out_msg)
        return False, out_msg

    passed = True
    failed_msg = ''

    # check pixel shift header keys for all files in obsdir
    for filename in tqdm(files, leave=False):
        # get snr
        snr = io.get_header_key(filename, SNR_KEY, dtype=float,
                                    required=False, default=None)
        # get dpr type
        dprtype = io.get_header_key(filename, 'DPRTYPE',
                                    required=False, default=None)
        # get object name
        objname = io.get_header_key(filename, 'DRSOBJN',
                                    required=False, default=None)
        # only check science observations
        if dprtype not in SCI_DPRTYPES:
            continue
        # fail
        if snr < SNR_LIMIT:
            # there is a shift
            passed = False
            failed_msg += (f'SNR[{SNR_KEY}] = {snr} (Should be > {SNR_LIMIT})'
                           f'\n\tDRSOBJN = {objname}'
                           f'\n\tFile: {filename}\n')

    if len(failed_msg) > 0:
        out_msg = (failed_msg)
    else:
        out_msg = (f'All observations SNR > {SNR_LIMIT}')
    if log:
        print(out_msg)
    return passed, out_msg


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
