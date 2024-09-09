#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-12 at 09:42

@author: cook
"""
from typing import Optional

import apero_checks
from apero_checks import override
from apero_checks.core import misc

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'apero_check_override.py'
# version, date, author
__VERSION__ = apero_checks.base.__VERSION__
__DATE__ = apero_checks.base.__DATE__
__AUTHOR__ = apero_checks.base.__AUTHOR__


# =============================================================================
# Define functions
# =============================================================================
def main(yaml_file: Optional[str] = None, obsdir: Optional[str] = None,
         test_name: Optional[str] = None, today: bool = False):
    # print splash
    apero_checks.splash('APERO Override checks')
    # get params updated for input yaml file
    all_params = apero_checks.load_params(yaml_file, obsdir, test_name, today)
    # loop around profiles
    for profile in all_params:
        # get profile params
        params = all_params[profile]
        # add profile name to parameters
        params['apero profile name'] = profile
        # find test
        allowed, test_type = override.find_override_test(params)
        # if not allowed
        if not allowed:
            msg = '\tCannot override test="{0}"'
            margs = [params['test_name']]
            misc.log_msg(msg.format(*margs), level='info')
            continue
        # run the tests
        test_results, overrides = override.override_tests(params,
                                                          test_type=test_type)
        # upload the tests
        apero_checks.upload_tests(params, test_results, test_type=test_type)
        # overrides must be saved
        apero_checks.store_overrides(params, overrides, test_type=test_type)
    # finish with an end message
    apero_checks.end_msg()


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------:
    # run main
    main()

# =============================================================================
# End of code
# =============================================================================
