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
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from astropy import units as uu
from astropy.time import Time, TimeDelta

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
def get_files_prefix(path: str, suffix: str,
                     last: Optional[Time] = None,
                     first: Optional[Time] = None
                     ) -> Tuple[List[str], List[str]]:
    """
    Get the prefixes for all files in a directory with a given suffix

    :param path: str, the path to the directory
    :param suffix: str, the suffix to search for


    :return: tuple, 1. list of file prefixes, 2. list of file paths
    """
    # storage of prefixes and paths
    file_prefixes, file_paths = [], []
    # rejected dirs
    rejected_dirs = []
    # walk around all subdirectories with root = path
    for root, dirs, files in os.walk(path, followlinks=True):
        # loop around files
        for filename in files:
            # get obsdir
            obsdir = str(root)
            # make sure it is a file
            if filename.endswith(suffix):
                # do not add files for observation directories that were
                #    previous rejected
                if obsdir in rejected_dirs:
                    continue
                # get the absolute path
                abspath = os.path.join(root, filename)
                # get the last modified time
                if first is not None and last is not None:
                    last_mod = os.path.getmtime(abspath)
                # not required
                else:
                    last_mod = 0
                # only keep files in the time range
                if first is not None and last_mod < first.unix:
                    # add to rejected (to skip check in future
                    rejected_dirs.append(obsdir)
                    continue
                elif last is not None and last_mod > last.unix:
                    # add to rejected (to skip check in future
                    rejected_dirs.append(obsdir)
                    continue
                # make sure it is a file
                file_prefixes.append(filename.split(suffix)[0])
                file_paths.append(abspath)
    # return the prefixes and paths
    return file_prefixes, file_paths


def get_last_file_time(path: str) -> Union[Time, None]:
    """
    Get the newest file in a directory

    :param path: str, the absolute path to the directory

    :return:
    """
    while True:
        # list all files in this path
        files = glob.glob(os.path.join(path, '*.fits'))
        # deal with no files
        if len(files) == 0:
            directories = glob.glob(os.path.join(os.path.dirname(path), '*'))
            directories.sort(key=os.path.getmtime)
            # deal with no directories
            if len(directories) == 0:
                return None
            # set the path to the last directory
            path = directories[-1]
        else:
            break

    # sort these by last modified time
    files.sort(key=os.path.getmtime)
    # get the last file
    last_file = files[-1]
    # return the last modified time
    return Time(os.path.getmtime(last_file), format='unix')


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
    apero_params = apero_functions.update_apero_profile(params)
    # define parameters we use here
    raw_directory = apero_params['DRS_DATA_RAW']
    tmp_directory = apero_params['DRS_DATA_WORKING']
    # -------------------------------------------------------------------------
    # get the last modified time of the last file in the raw directory +
    #   obs dir
    last_modified = get_last_file_time(os.path.join(raw_directory, obsdir))
    # deal with no last_modified --> raw directory is empty
    if last_modified is None:
        if log:
            msg = 'Raw directory is empty: {0}'.format(raw_directory)
            misc.log_msg(msg, level='warning')
        return False
    # -------------------------------------------------------------------------
    # get the first and last time (we check for the last 7 days only)
    first = last_modified - TimeDelta(7 * uu.day)
    last = last_modified
    # -------------------------------------------------------------------------
    # get all raw files prefix
    raw_prefixes, raw_files = get_files_prefix(raw_directory, '.fits',
                                               first=first, last=last)
    pp_prefixes, _ = get_files_prefix(tmp_directory, '_pp.fits')
    # -------------------------------------------------------------------------
    # loop around all files and find any raw files that are not in pp
    no_reduction = ~np.in1d(raw_prefixes, pp_prefixes)
    # if we have any files that are not in reduction we log an error
    if np.sum(no_reduction) > 0:
        if log:
            msg = 'Missing pp files in last 7 days: '
            for row in range(len(raw_prefixes)):
                if no_reduction[row]:
                    msg += '\n\t{0}'.format(raw_files[row])

            misc.log_msg(msg, level='warning')
        return False
    # -------------------------------------------------------------------------
    if log:
        msg = ('No missing pp files in the last 7 days.')
        misc.log_msg(msg, level='')
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
