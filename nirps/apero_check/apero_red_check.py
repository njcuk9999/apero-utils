#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
from typing import Optional

import apero_checks

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = apero_checks.base.__VERSION__
__DATE__ = apero_checks.base.__DATE__
__AUTHOR__ = apero_checks.base.__AUTHOR__


# -----------------------------------------------------------------------------


# =============================================================================
# Define functions
# =============================================================================
def main(yaml_file: Optional[str] = None, obsdir: Optional[str] = None,
         test_name: Optional[str] = None, today: bool = False,
         yest: bool = False, since: str = 'None', until: str = 'None'):
    # print splash
    apero_checks.splash('APERO Reduction checks')
    # get params updated for input yaml file
    all_params = apero_checks.load_params(yaml_file, obsdir, test_name, today,
                                          yest, since, until)
    # set up log results
    log_results = apero_checks.define_log_results()
    # loop around profiles
    for profile in all_params:
        # get profile params
        params = all_params[profile]
        # add profile name to parameters
        params['apero profile name'] = profile
        # if we do not have a test name then we run all tests and upload
        if params['test_name'] in [None, 'None']:
            # run the tests
            test_results = apero_checks.run_tests(params, log_results,
                                                  test_type='red')
            # upload the tests
            apero_checks.upload_tests(params, test_results, test_type='red')
        # otherwise we run a single test
        else:
            # run single test
            apero_checks.run_single_test(params, log_results, test_type='red')
        # log failures to file
        apero_checks.log_tests(log_results)
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
