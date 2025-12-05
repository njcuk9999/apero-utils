#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import smtplib
import shutil
from typing import Any, Dict, List, Optional, Tuple, Union

from apero_checks.core import misc


# =============================================================================
# Define variables
# =============================================================================
DISK_PATH = '/cosmos99'
# Cache results (no need to repeat)
CACHE = dict()

# =============================================================================
# Define functions
# =============================================================================
def test_email_server(params) -> Tuple[bool, str]:
    """
    Test the email server is working (used to send emails)

    :return:
    """
    # disk path
    email_server = params['system email server']
    if email_server == 'None':
        msg = 'Email server check disabled'
        print('"check.system email server = None"'
              '\n\t - No email server set, skipping test')
        return True, msg

    try:
        with smtplib.SMTP(email_server) as smtp:
            smtp.ehlo()  # Identify ourselves
            reason = ('Email server is reachable and responding')
            return True, reason
    except Exception as e:
        emsg = 'Failed to connect to email server\n\t{0}:{1}'
        eargs = [type(e), str(e)]
        return False, emsg.format(*eargs)


def test_disk_space(params) -> Tuple[bool, str]:
    """
    Test the disk space is sufficient for reductions
    :return:
    """
    # disk path
    disk_path = params['system disk path']
    if disk_path == 'None':
        msg = 'Disk usage check disabled'
        print('"check.system email server = None"'
              '\n\t - No disk path set, skipping test')
        return True, msg
    # get the disk stats
    total, used, free = shutil.disk_usage(disk_path)
    # calculate the disk usage
    usage = (used / total) * 100
    # fail if usage it over 90%
    if usage > 90:
        msg = 'Disk usage: {0} > 90% (Currently = {1:.2f}%)'
        margs = [disk_path, usage]
        return False, msg.format(*margs)
    else:
        msg = 'Disk usage: {0} okay (Currently = {1:.2f}%)'
        margs = [disk_path, usage]
        return True, msg.format(*margs)


def test(params: Dict[str, Any], obsdir: str, log=False) -> Tuple[bool, str]:
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
    # allow modification of the cache
    global CACHE
    # don't currently use params or obsdir
    _ = params, obsdir
    # all tests pass by default
    passed = True
    # define sub-tests
    test_functions = dict()
    test_functions['EMAIL SERVER'] = test_email_server
    test_functions['DISK_USAGE'] = test_disk_space

    # -------------------------------------------------------------------------
    # run all tests
    # -------------------------------------------------------------------------
    all_out_msg = ''
    # loop around tests
    for it, test_name in enumerate(list(test_functions.keys())):
        # if we have a cached result use it (we don't need to run these tests
        # multiple times)
        if test_name in CACHE:
            pass_it, reason_it = CACHE[test_name]
        else:
            # run this test
            pass_it, reason_it = test_functions[test_name](params)
            # add to cache
            CACHE[test_name] = (pass_it, reason_it)
        # print progress
        if log:
            misc.log_msg(f'Running {test_name} test [{it+1} of '
                         f'{len(test_functions)}]', '')
        # deal with failure
        out_msg = f'\t{reason_it}'

        if log and pass_it:
            misc.log_msg(out_msg, level='warning')
        elif log:
            misc.log_msg(out_msg, level='')

        all_out_msg += '\n' + out_msg
        # update passed
        passed &= pass_it
    # -------------------------------------------------------------------------
    return passed, all_out_msg


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
