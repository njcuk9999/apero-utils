#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:36

@author: cook
"""
# only import from this directory
from apero_checks.raw_tests import blank_test
from apero_checks.raw_tests import calib_test
from apero_checks.raw_tests import eng_test
from apero_checks.raw_tests import obsdir_test
from apero_checks.raw_tests import qual_test
from apero_checks.raw_tests import astrom_test
from apero_checks.raw_tests import prev_sci_test

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
test_dict['HAS_OBSDIR'] = obsdir_test.test

# eng test - Check the consistency of headers with what we expect from
#            enginnering. Consistency of FP temperature with set point,
#            enclosure, pumps, valves
test_dict['ENG_TEST'] = eng_test.test

# eng test - Are the expected calibrations present in the night folder
#            add eventual extra calibrations to the 'must have' objects
test_dict['CALIB_TEST'] = calib_test.test

# science quality test - Some basic quality checks for science files. Currently:
#                        saturation, flux
test_dict['SCI_QUAL'] = qual_test.sci_qual_test

# calib quality test - Some basic quality checks for calibration files. Currently:
#                      saturation, flux
test_dict['CALIB_QUAL'] = qual_test.calib_qual_test

# astrometric object test - If we have any objects missing from the astrometric
#                           database this should return False
test_dict['ASTROM_TEST'] = astrom_test.test


# previous science data test - this tests whether the previous 3 nights had
#                              science data
test_dict['PREV_SCI'] = prev_sci_test.test

# test 1 - explanation

# test 2 - explanation

# test 3 - explanation



# =============================================================================
# If and only if you want the user to be able to override test
# =============================================================================
# dictionary to store all tests
override_list = []

# PREV_SCI: There could be nights which do not have science data but should
#           not be flagged as bad
override_list.append('PREV_SCI')


# =============================================================================
# End of code
# =============================================================================
