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
from apero_checks.red_tests import manual_trigger_tests
from apero_checks.red_tests import pixel_shifts_test
from apero_checks.red_tests import excess_modal
from apero_checks.red_tests import prev_reduc_test

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
test_dict['HAS_OBSDIR'] = obsdir_test.test

# APERO calib test - test the APERO requirements on calibration
test_dict['APERO_CALIB'] = calib_test.test

# Test that manual trigger was started
test_dict['MANUAL_START'] = manual_trigger_tests.test_manual_trigger_start

# Test that manual trigger ended
test_dict['MANUAL_END'] = manual_trigger_tests.test_manual_trigger_end

# Test that APERO started
test_dict['APERO_START'] = manual_trigger_tests.test_apero_start

# Test that APERO ended
test_dict['APERO_END'] = manual_trigger_tests.test_apero_end

# Test that LBL started
test_dict['ARI_START'] = manual_trigger_tests.test_ari_start

# Test that LBL ended
test_dict['ARI_END'] = manual_trigger_tests.test_ari_end

# Pixel shifts tests - check for pixel shifts in pp files
test_dict['PIXEL_SHIFTS'] = pixel_shifts_test.test

# Test excess modal noise in telluric stars
test_dict['EXCESS_MODAL'] = excess_modal.test

# Test that every raw file has a preprocessed file
test_dict['PREV_REDUC'] = prev_reduc_test.test

# test 1 - explanation

# test 2 - explanation

# test 3 - explanation

# =============================================================================
# If and only if you want the user to be able to override test
# =============================================================================
# dictionary to store all tests
override_list = []


# =============================================================================
# End of code
# =============================================================================
