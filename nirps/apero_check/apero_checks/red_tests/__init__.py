#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:36

@author: cook
"""
# only import from this directory
from apero_checks.red_tests import blank_test
from apero_checks.red_tests import obsdir_test
from apero_checks.red_tests import calib_test

# =============================================================================
# Append your test here
# =============================================================================
# dictionary to store all tests
test_dict = dict()

# form of test_dict entry is as follows:
#    test_dict['KEY'] = test_module.function
#
#    where function must be in form
#    function(params: Dict[str, Any], obsdir: str, log=False) -> bool:
#
# Note the key is the column name that is used
#    Please do not have punctuation (other than underscore) and no spaces

# blank test - this tests whether the apero red tests ran
test_dict['BLANK'] = blank_test.test

# plan of tests
# - apero calibration requirements: use functionality from apero_precheck
# - apero astrometrics: object names (how to do per night?)
# - apero reduction ran on night
# - apero log errors
# - apero log qc passes
# - lbl ran on night
# - lbl errors?

# Manual trigger has run
test_dict['OBSDIR'] = obsdir_test.test

# APERO calib test - test the APERO requirements on calibration
test_dict['APERO_CALIB'] = calib_test.test

# test 1 - explanation

# test 2 - explanation

# test 3 - explanation


# =============================================================================
# End of code
# =============================================================================
