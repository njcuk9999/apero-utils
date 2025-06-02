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
import numpy as np

from apero_checks.core import apero_functions
from apero_checks.core import misc
from apero_checks.core import io


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)

# science DPRTYPES
SCI_DPRTYPES = ['OBJ_DARK', 'OBJ_FP', 'OBJ_SKY']
# set ccf out file
CCF_OUT_FILE = "CCF_RV"
# Threshold for flagging bad files (in units of sigma)
BAD_NSIG = 10.0


# =============================================================================
# Define functions
# =============================================================================
def sigma(x):
    """Calculate the robust standard deviation of a dataset."""
    # Median Absolute Deviation (MAD) scaled to approximate std
    # deviation for normal dist.
    return np.nanmedian(np.abs(x - np.nanmedian(x))) * 1.4826


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
        out_msg = ('red directory {} does not exist'.format(obsdir))
        if log:
            print(out_msg)
        return False, out_msg
    # -------------------------------------------------------------------------
    # import apero here
    from apero.core.core import drs_database
    # get file index database
    findexdb = drs_database.FileIndexDatabase(apero_params)
    findexdb.load_db()
    # -------------------------------------------------------------------------
    # Get a list of targets for this night
    # -------------------------------------------------------------------------
    # look only at raw data
    condition = 'BLOCK_KIND="raw" AND OBS_DIR="{obsdir}"'
    # look only at science dprtypes
    subconds = []
    for dprtype in SCI_DPRTYPES:
        subconds.append(f'KW_DPRTYPE="{dprtype}"')
    condition += ' AND ({0})'.format(' OR '.join(subconds))
    # run query
    sci_objnames = findexdb.get_unique('KW_OBJNAME', condition=condition)
    # deal with no science objects this night (return True)
    if len(sci_objnames) == 0:
        out_msg = 'No SCI objects found for {0}'.format(obsdir_path)
        return True, out_msg
    elif log:
        print(' Found {0} SCI objects'.format(len(sci_objnames)))

    # -------------------------------------------------------------------------
    # loop around each object
    # -------------------------------------------------------------------------
    # storage for return
    out_msg = ''
    passed = True
    # loop around objects
    for objname in tqdm(sci_objnames):
        # get ccf condition
        condition = f'BLOCK_KIND="red" AND KW_OUTPUT="{CCF_OUT_FILE}"'
        condition += f' AND KW_OBJNAME="{objname}"'

        # query database for filenames
        filenames = findexdb.get_entries('ABSPATH', condition=condition)

        # get rv_obj and ccf_mfwhm from headers
        rv_obj, ccf_mfwhm, valid_files = [], [], []
        for filename in filenames:
            if not os.path.exists(filename):
                continue
            # read header
            hdr = fits.getheader(filename)
            # get values
            rv_obj.append(float(hdr['RV_OBJ']))
            ccf_mfwhm.append(float(hdr['CCFMFWHM']))
        # convert to numpy arrays
        rv_obj = np.array(rv_obj)
        ccf_mfwhm = np.array(ccf_mfwhm)
        valid_files = np.array(valid_files)

        # Compute median and robust std deviation for each parameter
        median_rv = np.nanmedian(rv_obj)
        median_fwhm = np.nanmedian(ccf_mfwhm)
        sig_rv = sigma(rv_obj)
        sig_fwhm = sigma(ccf_mfwhm)

        # Calculate the "distance" in sigma units for each file from the median
        part1 = ((rv_obj - median_rv) / sig_rv)
        part2 = ((ccf_mfwhm - median_fwhm) / sig_fwhm)
        nsig = np.sqrt(part1** 2 +  part2**2)

        # Identify files that are farther than bad_nsig from the
        # median (outliers)
        bad = nsig > BAD_NSIG

        bad_files = valid_files[bad]
        # if we have bad files state it here
        if bad < 0:
            passed = False
            # add to out message
            out_msg += f'FAILED: {objname} outlier CCFs identified'
            # loop around bad files
            for it, bad_file in enumerate(bad_files):
                out_msg += f'\t\t{it + 1}: {bad_file}'
            out_msg += '\nn'
        # otherwise we don't have bad files
        else:
            out_msg += f'PASSED: {objname} had no bad outliers\n\n'

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
