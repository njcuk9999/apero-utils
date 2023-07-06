#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:36

@author: cook
"""
# only import from this directory
from apero_raw_tests.tests import blank_test
from apero_raw_tests.tests import obsdir_test
from apero_raw_tests.tests import eng_test
from apero_raw_tests.tests import calib_test

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

# blank test - this tests whether the apero raw tests ran
test_dict['BLANK'] = blank_test.test

# obs dir test - this tests whether the obsdir given exists on disk and has
#                fits files in its directory
test_dict['OBSDIR'] = obsdir_test.test

# eng test - Check the consistency of headers with what we expect from
#            enginnering. Consistency of FP temperature with set point,
#            enclosure, pumps, valves
test_dict['ENG_TEST'] = eng_test.test

# eng test - Are the expected calibrations present in the night folder
#            add eventual extra calibrations to the 'must have' objects
test_dict['CALIB_TEST'] = calib_test.test

# test 1 - explanation

# test 2 - explanation

# test 3 - explanation



# =============================================================================
# End of code
# =============================================================================
