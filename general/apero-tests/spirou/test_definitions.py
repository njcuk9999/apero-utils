"""
Definition of SPIRou recipe tests following APERO framework

@author: vandalt
"""
import sys  # Temporary

from apero.core import constants
from apero.core.instruments.spirou.recipe_definitions import recipes


# Temporary to import full paths
sys.path.insert(0, '..')
from drs_test import DrsTest  # Eventually from apero-test.drs_test

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = "apero-test.spirou.test_definitions.py"
__INSTRUMENT__ = "SPIROU"
# Get constants
Constants = constants.load(__INSTRUMENT__)
# Get version and author
__version__ = Constants["DRS_VERSION"]
__author__ = Constants["AUTHORS"]
__date__ = Constants["DRS_DATE"]
__release__ = Constants["DRS_RELEASE"]

RECIPE_DICT = dict(zip(list(map(lambda x: x.name, recipes)), recipes))

# =============================================================================
# Define tests
# =============================================================================
tests = []

# -----------------------------------------------------------------------------
# Badpixel Test
# -----------------------------------------------------------------------------
cal_badpix = DrsTest(__INSTRUMENT__, RECIPE_DICT["cal_badpix_spirou.py"])
cal_badpix.name = "Bad Pixel Correction Recipe Test #1"
# TODO: Add other attributes here
tests.append(cal_badpix)

# TODO: Add other tests
