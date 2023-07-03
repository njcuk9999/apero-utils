#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
from typing import Optional

import apero_raw_tests

# =============================================================================
# Define variables
# =============================================================================
__VERSION__ = apero_raw_tests.base.__VERSION__
__AUTHOR__ = apero_raw_tests.base.__AUTHOR__
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def main(yaml_file: Optional[str] = None, obsdir: Optional[str] = None,
         test_name: Optional[str] = None):
    # get params updated for input yaml file
    params = apero_raw_tests.load_params(yaml_file, obsdir)
    # if we do not have a test name then we run all tests and upload
    if params['test_name'] not in [None, 'None']:
        # run the tests
        test_results = apero_raw_tests.run_tests(params)
        # upload the tests
        apero_raw_tests.upload_tests(params, test_results)
    # otherwise we run a single test
    else:
        # run single test
        apero_raw_tests.run_single_test(params, test_name)


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
