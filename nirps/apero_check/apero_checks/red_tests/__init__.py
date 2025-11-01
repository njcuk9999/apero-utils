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
from apero_checks.red_tests import low_snr_test
from apero_checks.red_tests import ccf_test
from apero_checks.red_tests import prev_reduc_test


# =============================================================================
# Append your test here
# =============================================================================
# dictionary to store all tests
test_dict = dict()
test_dep = dict()

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
test_dep['BLANK'] = []

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
test_dep['HAS_OBSDIR'] = ['BLANK']

# APERO calib test - test the APERO requirements on calibration
test_dict['APERO_CALIB'] = calib_test.test
test_dep['APERO_CALIB'] = ['BLANK', 'HAS_OBSDIR']

# Test that manual trigger was started
test_dict['MANUAL_START'] = manual_trigger_tests.test_manual_trigger_start
test_dep['MANUAL_START'] = ['BLANK', 'HAS_OBSDIR']

# Test that manual trigger ended
test_dict['MANUAL_END'] = manual_trigger_tests.test_manual_trigger_end
test_dep['MANUAL_END'] = ['BLANK', 'HAS_OBSDIR', 'MANUAL_START']

# Test that APERO started
test_dict['APERO_START'] = manual_trigger_tests.test_apero_start
test_dep['APERO_START'] = ['BLANK', 'HAS_OBSDIR', 'MANUAL_START']

# Test that APERO ended
test_dict['APERO_END'] = manual_trigger_tests.test_apero_end
test_dep['APERO_END'] = ['BLANK', 'HAS_OBSDIR', 'MANUAL_START']

# Test that LBL started
test_dict['ARI_START'] = manual_trigger_tests.test_ari_start
test_dep['ARI_START'] = ['BLANK', 'HAS_OBSDIR', 'MANUAL_START']

# Test that LBL ended
test_dict['ARI_END'] = manual_trigger_tests.test_ari_end
test_dep['ARI_END'] = ['BLANK', 'HAS_OBSDIR', 'MANUAL_START', 'ARI_START']

# Pixel shifts tests - check for pixel shifts in pp files
test_dict['PIXEL_SHIFTS'] = pixel_shifts_test.test
test_dep['PIXEL_SHIFTS'] = ['BLANK', 'HAS_OBSDIR', 'APERO_START', 'APERO_END']

# Test excess modal noise in telluric stars
test_dict['EXCESS_MODAL'] = excess_modal.test
test_dep['EXCESS_MODAL'] = ['BLANK', 'HAS_OBSDIR', 'APERO_START', 'APERO_END']

# Test for low SNR in science targets
test_dict['LOW_SNR'] = low_snr_test.test
test_dep['LOW_SNR'] = ['BLANK', 'HAS_OBSDIR', 'APERO_START', 'APERO_END']

# Test for bad CCF in science targets
test_dict['BAD_CCF'] = ccf_test.test
test_dep['BAD_CCF'] = ['BLANK', 'HAS_OBSDIR', 'APERO_START', 'APERO_END']

# Test that every raw file has a preprocessed file
test_dict['PREV_REDUC'] = prev_reduc_test.test
test_dep['PREV_REDUC'] = ['BLANK', 'HAS_OBSDIR', 'APERO_START', 'APERO_END']

# test 1 - explanation

# test 2 - explanation

# test 3 - explanation

# =============================================================================
# If and only if you want the user to be able to override test
# =============================================================================
# dictionary to store all tests
override_list = []

# LOW_SNR: There could be objects for which we expect a low SNR
override_list.append('LOW_SNR')

# =============================================================================
# End of code
# =============================================================================
