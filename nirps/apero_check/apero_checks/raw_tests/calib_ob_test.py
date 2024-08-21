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
from typing import Any, Dict

from astropy.io import fits
from tqdm import tqdm

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
    Test with observation directory exists

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # Here we add OB NAMES (sub strings) that are required for a night
    #  note these is MJD sensitive so we require a range of MJDs
    must_have_obnames = dict()
    # Afternoon calibrations started around 2023-06-01 - only check after
    must_have_obnames['Afternoon-Calibrations'] = (60096, 99999)
    # Daily calibration started around 2022-11-12 - only check after
    must_have_obnames['Daily-Calibrations-Full-Set'] = (59895, 99999)
    # This tells use what the header keys are
    header_keys = dict()
    header_keys['DPR_TYPE'] = 'HIERARCH ESO DPR TYPE'
    header_keys['OBS_NAME'] = 'HIERARCH ESO OBS NAME'
    header_keys['MJD-OBS'] = 'MJD-OBS'
    # define parameters we use here
    raw_directory = params['raw dir']
    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------
    # get the observation directory
    obsdir_path = os.path.join(raw_directory, obsdir)
    # deal with no observation directory
    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {0} does not exist - TEST FAILED')
        return False
    # get the files for this observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))
    # deal with no files
    if len(files) == 0:
        if log:
            # pass a True if no file is found on that night
            print('No files found for night {}'.format(obsdir))
        return False
    # -------------------------------------------------------------------------
    # create table to store keywords
    obs_counts = dict()
    # flag obs out of range
    obs_out_of_range = dict()

    # fill table looping through files
    for ifile in tqdm(range(len(files))):
        # get the header
        hdr = fits.getheader(files[ifile])
        # get dpr type and obs name
        obs_name = str(hdr.get(header_keys['OBS_NAME'], None))
        mjd_date = hdr[header_keys['MJD-OBS']]
        # check must have obnames
        for obname in must_have_obnames:
            # get the conditions to add obname
            in_range = mjd_date >= must_have_obnames[obname][0]
            in_range &= mjd_date <= must_have_obnames[obname][1]
            # if we are out of range we don't test this obname
            if not in_range:
                obs_out_of_range[obname] = True
            elif obname in obs_name:
                # count obs name instances
                if obs_name in obs_counts:
                    obs_counts[obname] += 1
                else:
                    obs_counts[obname] = 1
        # close the header
        del hdr
    # -------------------------------------------------------------------------
    # set the failed messages
    failed_log = []
    passed_log = []
    # -------------------------------------------------------------------------
    # Check that we have morning and afternoon calibrations
    # -------------------------------------------------------------------------
    # loop around obnames
    for obs_name in must_have_obnames:
        # do not test if out of range
        if obs_name in obs_out_of_range:
            continue
        # deal with pass/failure
        if obs_name in obs_counts:
            passed_log.append('\t"{0}" OBNAME found.'.format(obs_name))
        else:
            failed_log.append('\t"{0}" OBNAME not found'.format(obs_name))
    # -------------------------------------------------------------------------
    # construct a string for printing output
    passed_log = '\n'.join(passed_log)
    if len(passed_log) == 0:
        passed_log = '\tNone'
    failed_log = '\n'.join(failed_log)

    if len(failed_log) == 0:
        failed_log = '\tNone'
        passed_logical = True
    else:
        passed_logical = False
    # -------------------------------------------------------------------------
    if log:
        print('\n')
        print('*' * 50)
        print('QC for night {}'.format(obsdir))
        print('Passed QC:')
        print(passed_log)
        print('Failed QC:')
        print(failed_log)
        print('*' * 50)
        print('\n')
    # -------------------------------------------------------------------------
    return passed_logical


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-07-01'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
