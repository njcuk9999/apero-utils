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


def sat_test(files, log=False) -> bool:
    """
    Saturation test - check for abnormal saturation in raw files.
    Uses the n_read (ext=3) extension to compute saturation fraction.
    With log on, print failed files and the failure reason.
    Called by the main calib/sci_qual_test functions

    All files passed = True
    Any file failed = False

    :param files: a list of files of a single DPR type, obtained in the
                  main test function
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """

    failed_log = ''

    for filename in tqdm(files, leave=False):

        dpr_type = io.get_header_key(filename, 'HIERARCH ESO DPR TYPE')

        basename = os.path.basename(filename)

        try:
            constraints_type = constraints[dpr_type]
        except KeyError:
            continue
        nread = np.array(fits.getdata(filename, ext=3))

        # saturation fraction
        fsat = np.mean(nread != np.max(nread))

        if fsat > constraints_type[2]:
            args = [basename, constraints_type[2], fsat]
            failed_log += ('{}: fsat is too high, limit at {}, '
                           'value at {:.2f}\n').format(*args)

    passed = True
    if len(failed_log) > 0:
        passed = False

    if log:
        print('Saturation check:')
        if len(failed_log) > 0:
            print('Issues were found in the following files:')
            print(failed_log)
        else:
            print('All passed.\n')

    return passed


def flux_test(files, log=False) -> bool:
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

    for filename in tqdm(files, leave=False):

        dpr_type = io.get_header_key(filename, 'HIERARCH ESO DPR TYPE')

        basename = os.path.basename(filename)

        try:
            constraints_type = constraints[dpr_type]
        except KeyError:
            continue
        im = np.array(fits.getdata(filename, ext=1))

        # 99th percentile
        p99 = np.nanpercentile(im, 99)

        if p99 < constraints_type[0]:
            args = [basename, constraints_type[0], p99]
            failed_log += ('{}: p99 is too low, limit at {}, '
                           'value at {:.2f}\n').format(*args)
        if p99 > constraints_type[1]:
            args = [basename, constraints_type[1], p99]
            failed_log += ('{}: p99 is too high, limit at {}, '
                           'value at {:.2f}\n').format(*args)

    passed = True
    if len(failed_log) > 0:
        passed = False

    if log:
        print('Flux check:')
        if len(failed_log) > 0:
            print('Issues were found in the following files:')
            print(failed_log)
        else:
            print('All passed.\n')

    return passed


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

    # select FP, FLAT and TELLURIC files based on DPR types
    fp_files = []
    flat_files = []
    telluric_files = []
    for filename in tqdm(files, leave=False):

        dpr_type = io.get_header_key(filename, 'HIERARCH ESO DPR TYPE')

        if dpr_type in dpr_types['FP']:
            fp_files.append(filename)
        if dpr_type in dpr_types['FLAT']:
            flat_files.append(filename)
        if dpr_type in dpr_types['TELLURIC']:
            telluric_files.append(filename)

    # run the tests for each DPR type
    outputs = []

    dpr_tuple = [('FP', fp_files), ('FLAT', flat_files),
                 ('TELLURIC', telluric_files)]

    for (dpr_type, dpr_files) in dpr_tuple:
        if log:
            print('{} FILES:\n'.format(dpr_type))
        outputs.append(sat_test(dpr_files, log=log))
        outputs.append(flux_test(dpr_files, log=log))

    # False if any of the tests for any file failed
    passed = True
    for output in outputs:
        passed = passed and output

    return passed


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

    # select only the science files
    science_files = []
    for filename in tqdm(files, leave=False):

        dpr_type = io.get_header_key(filename, 'HIERARCH ESO DPR TYPE')

        if dpr_type in dpr_types['SCIENCE']:
            science_files.append(filename)

    # run the tests on the science files
    if log:
        print('SCIENCE FILES:\n')
    sat_passed = sat_test(science_files, log=log)
    flux_passed = flux_test(science_files, log=log)

    passed = sat_passed and flux_passed

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
    calib_qual_test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
