"""
Defintion of SPIRou recipes following APER framework
"""
from drs_test import DrsTest

# =============================================================================
# Define tests
# =============================================================================
tests = []

# -----------------------------------------------------------------------------
# Badpixel Test
# -----------------------------------------------------------------------------
cal_badpix = DrsTest()
cal_badpix.long_name = 'Bad Pixel Correction Recipe Test #1'  # Old self.name
# TODO: Add other attributes here
tests.append(cal_badpix)

# TODO: Add other tests
