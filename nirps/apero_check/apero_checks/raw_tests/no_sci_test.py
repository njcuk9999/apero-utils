#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import glob
import os
from typing import Any, Dict, Tuple

from tqdm import tqdm

from apero_checks.core import io

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)

# dprtype header key
DPRTYPE_KEY = 'HIERARCH ESO DPR TYPE'

# dprtype key that defines science files
SCIENCE_DPRTYPES = ['OBJECT,SKY', 'OBJECT,FP']

# We define this manually as we need to check HE and HA
PATHS_TO_RAW = ['/nirps_raw/nirps/raw-data/nirps_ha',
                '/nirps_raw/nirps/raw-data/nirps_he']

# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> Tuple[bool, str]:
    """
    Blank test - this tests whether test was run (should always return True)
    All other tests should return True or False, and only print messages if
    log is True.

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # we don't use params
    _ = params
    # define parameters we use here
    files = []
    # loop around all raw directories
    for raw_directory in PATHS_TO_RAW:
        # get a list of all files
        obsdir_path = os.path.join(raw_directory, obsdir)
        # skip if path doesn't exist
        if not os.path.exists(obsdir_path):
            continue
        # get all files for this night
        files += glob.glob(os.path.join(obsdir_path, '*.fits'))
    # no files is bad
    if len(files) == 0:
        # pass a True if no file is found on that night
        out_msg = ('No files found for night {}'.format(obsdir))
        for raw_directory in PATHS_TO_RAW:
            # get a list of all files
            obsdir_path = os.path.join(raw_directory, obsdir)
            # print paths checked
            out_msg += ('\n\tChecked: {0}'.format(obsdir_path))
        if log:
            print(out_msg)
        # no files mean we fail
        return False, out_msg
    # -------------------------------------------------------------------------
    # storage for counter
    science_files = dict()
    other_files = dict()
    # -------------------------------------------------------------------------
    # loop around files
    for ifile in tqdm(range(len(files))):
        # get file for this iteration
        filename = files[ifile]
        # get dprtype
        dprtype = io.get_header_key(filename, DPRTYPE_KEY,
                                      required=False, default=None)
        # check if we have a science file
        if dprtype in SCIENCE_DPRTYPES:
            # increment counter if we already have this file
            if dprtype in science_files:
                science_files[dprtype] += 1
            # else set counter to 1
            else:
                science_files[dprtype] = 1
        else:
            # increment counter if we already have this file
            if dprtype in other_files:
                other_files[dprtype] += 1
            # else set counter to 1
            else:
                other_files[dprtype] = 1
    # -------------------------------------------------------------------------

    # no science files might be bad
    if len(science_files) == 0:
        out_msg = ('No science files found for night {}'.format(obsdir))
        out_msg += ('\nOther files found for night {}'.format(obsdir))
        for key in other_files:
            out_msg += ('\n\t{0}: {1}'.format(key, other_files[key]))
        out_msg += ('\n\nNote: If override used it will not display here. ')
        if log:
            print(out_msg)

        return False, out_msg
    # -------------------------------------------------------------------------
    # log number of science files found for each type
    out_msg = ('Science files found for night {0}'.format(obsdir))
    for key in science_files:
        out_msg += ('\n\t{0}: {1}'.format(key, science_files[key]))
    out_msg += ('\n\tOther files found for night {}'.format(obsdir))
    for key in other_files:
        out_msg += ('\n\t{0}: {1}'.format(key, other_files[key]))
    if log:
        print(out_msg)
    return True, out_msg


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
