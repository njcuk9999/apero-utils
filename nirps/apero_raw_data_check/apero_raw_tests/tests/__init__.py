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

# =============================================================================
# Append your test here
# =============================================================================
test_dict = dict()
# blank test - this tests whether the apero raw tests ran
test_dict['BLANK'] = blank_test.test
# test 1
test_dict['OBSDIR'] = obsdir_test.test
# test 2
test_dict['ENG_TEST'] = eng_test.test
# test 3



# =============================================================================
# End of code
# =============================================================================
