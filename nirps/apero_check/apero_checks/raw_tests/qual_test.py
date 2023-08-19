#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict, List, Tuple
import os
import glob
from tqdm import tqdm

from astropy.io import fits
import numpy as np

from apero_checks.core import misc
from apero_checks.core import io

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================

# list of DPR TYPE for each file category
dpr_types = {
    'FP': ['WAVE,FP,FP', 'WAVE,UN1,UN1', 'WAVE,FP,UN1', 'WAVE,UN1,FP',
           'CONTAM,DARK,FP'],
    'FLAT': ['FLAT,LAMP,DARK', 'FLAT,DARK,LAMP', 'FLAT,LAMP,LAMP',
             'ORDERDEF,DARK,LAMP', 'ORDERDEF,LAMP,DARK'],
    'TELLURIC': ['TELLURIC,SKY'],
    'SCIENCE': ['OBJECT,SKY', 'OBJECT,FP']
}

# define normal constraints for minimum 99th percentile, maximum 99th percentile
# and max fraction of saturated pixels
# for each DPR TYPE
constraints = dict()
constraints['OBJECT,SKY'] = 0, 10000, 0.0004
constraints['OBJECT,FP'] = 0, 10000, 0.0004
constraints['ORDERDEF,DARK,LAMP'] = 100, 1000, 0.0004
constraints['ORDERDEF,LAMP,DARK'] = 100, 1000, 0.0004
constraints['TELLURIC,SKY'] = 0, 10000, 0.0004
constraints['FLAT,LAMP,DARK'] = 100, 1000, 0.0004
constraints['FLAT,DARK,LAMP'] = 100, 1000, 0.0004
constraints['FLAT,LAMP,LAMP'] = 100, 1000, 0.0004
constraints['CONTAM,DARK,FP'] = 10, 1000, 0.0004
constraints['WAVE,FP,FP'] = 10, 1000, 0.0004
constraints['WAVE,UN1,UN1'] = 5, 100, 0.0025
constraints['WAVE,FP,UN1'] = 5, 100, 0.0012
constraints['WAVE,UN1,FP'] = 5, 100, 0.0023


def sat_test(filename, nread, dprtype) -> Tuple[bool, str]:
    """
    Saturation test - check for abnormal saturation in raw files.
    Uses the n_read (ext=3) extension to compute saturation fraction.
    With log on, print failed files and the failure reason.
    Called by the main calib/sci_qual_test functions

    All files passed = True
    Any file failed = False

    :param files: a list of files of a single DPR type, obtained in the
                  main test function

    :return: bool, True if passed, False otherwise
    """
    basename = os.path.basename(filename)
    try:
        constraints_type = constraints[dprtype]
    except KeyError:
        return True, ''
    # saturation fraction
    fsat = np.mean(nread != np.max(nread))

    if fsat > constraints_type[2]:
        args = [basename, constraints_type[2], fsat]
        failed_log = ('Saturation check: {}: fsat is too high, limit at {}, '
                       'value at {:.2f}\n').format(*args)
        return False, failed_log

    return True, ''


def flux_test(filename, image, dprtype) -> Tuple[bool, str]:
    """
    Flux test - check for abnormal flux level in raw files.
    Uses the image (ext=1) extension to compute 99th percentile and compare
    with expected limits.
    With log on, print failed files and the failure reason.
    Called by the main calib/sci_qual_test functions.

    All files passed = True
    Any file failed = False

    :param files: a list of files of a single DPR type, obtained in the
                  main test function
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """

    failed_log = ''

    basename = os.path.basename(filename)

    try:
        constraints_type = constraints[dprtype]
    except KeyError:
        return True, ''

    # 99th percentile
    p99 = np.nanpercentile(image, 99)

    if p99 < constraints_type[0]:
        args = [basename, constraints_type[0], p99]
        failed_log += ('Flux check: {}: p99 is too low, limit at {}, '
                       'value at {:.2f}\n').format(*args)
    if p99 > constraints_type[1]:
        args = [basename, constraints_type[1], p99]
        failed_log += ('Flux check: {}: p99 is too high, limit at {}, '
                       'value at {:.2f}\n').format(*args)

    if len(failed_log) > 0:
        return False, failed_log

    return True, ''



def qual_test(params: Dict[str, Any], obsdir: str, dprgroups: List[str],
              log: bool = False) -> bool:
    # get path to observation directory
    raw_dir = params['raw dir']
    obsdir_path = os.path.join(raw_dir, obsdir)

    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    # check if directory exists
    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {} does not exist'.format(obsdir))
        return False

    # get a list of all the files in observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))

    # check if there are files in the directory
    if len(files) == 0:
        if log:
            print('No files in directory {}'.format(obsdir))
        return False
    # -------------------------------------------------------------------------
    # filter files to only look at those with the dprtypes we want
    # -------------------------------------------------------------------------
    # get all valid dprtypes
    valid_dprtypes = []
    for dprtype in dprgroups:
        valid_dprtypes += dpr_types[dprtype]
    # storage of valid files
    valid_files = []
    # loop around files and get valid files
    for filename in files:
        # get the dprtype
        dpr_type = fits.getheader(filename)['HIERARCH ESO DPR TYPE']
        # if dprtype is valid add to valid files
        if dpr_type in valid_dprtypes:
            valid_files.append(filename)
    # -------------------------------------------------------------------------
    # storage for all failed outputs
    failed_outputs = dict()
    failed_count = 0
    # make storage
    for dprtype in dprgroups:
        # each group is a dictionary of files
        failed_outputs[dprtype] = dict()
    # loop around all files
    for filename in tqdm(files, leave=False):
        # lets only open the file once
        with fits.open(filename) as hdul:
            image = hdul[1].data
            nread = hdul[3].data
            hdr = hdul[0].header
            # get dprtype
            dprtype = hdr['HIERARCH ESO DPR TYPE']
            for test_dprtype in dprgroups:
                if dprtype in dpr_types[test_dprtype]:
                    passed1, fail_msg1 = sat_test(filename, nread, dprtype)
                    passed2, fail_msg2 = flux_test(filename, image, dprtype)
                    # combine tests
                    passed = passed1 and passed2
                    fail_msgs = []
                    if passed1:
                        fail_msgs.append('\t' + fail_msg1)
                    if passed2:
                        fail_msgs.append('\t' + fail_msg2)

                    if not passed:
                        failed_outputs[test_dprtype][filename] = fail_msgs
                        failed_count += 1
                    break
    # -------------------------------------------------------------------------
    # False if any of the tests for any file failed
    if failed_count > 0:
        return False
    else:
        return True


def calib_qual_test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Calibration quality tests - multiple basic quality checks for
    calibration files.
    Checks 3 types of raws/DPR types: FP, FLAT, and TELLURIC
    With log on, print a list of files that failed each test by category/type

    Tests:
    - Saturation
    - Flux

    All files passed = True
    Any file failed = False

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # define the calibration dprtype groups
    dprgroups = ['FP', 'FLAT', 'TELLURIC']
    # run the quality tests
    return qual_test(params, obsdir, dprgroups, log=log)



def sci_qual_test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Calibration quality tests - multiple basic quality checks for science
    files.
    Checks 1 type of raws/DPR types: SCIENCE
    With log on, print a list of files that failed each test

    Tests:
    - Saturation
    - Flux

    All files passed = True
    Any file failed = False

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # define the calibration dprtype groups
    dprgroups = ['SCIENCE']
    # run the quality tests
    return qual_test(params, obsdir, dprgroups, log=log)


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
    calib_qual_test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
