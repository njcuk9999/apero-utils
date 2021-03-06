"""
Definition of SPIRou recipe tests following APERO framework

@author: vandalt
"""
import sys  # Temporary

from apero.core import constants
from apero.core.instruments.spirou.recipe_definitions import recipes

# HACK: Temporary to import full paths
sys.path.insert(0, '..')
from drs_test import DrsTest  # Eventually from apero-tests.drs_test

# Something similar?
# from .preprocessing import PPTest
# from .darkmaster import DarkMTest
from tests.badpixel import BadPixTest

# from .localisation import LocTest
# from .shapemaster import ShapeMTest
# from .shape import ShapeTest
# from .flat import FlatTest
# from .thermal import ThermalTest
# from .leakmaster import LeakMTest
# from .wavelengthmaster import WaveMTest
# from .wavelength import WaveTest
import test_utils as tu

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = "apero-tests.spirou.test_definitions.py"
__INSTRUMENT__ = "SPIROU"
# Get constants
Constants = constants.load(__INSTRUMENT__)
# Get version and author
__version__ = Constants["DRS_VERSION"]
__author__ = Constants["AUTHORS"]
__date__ = Constants["DRS_DATE"]
__release__ = Constants["DRS_RELEASE"]

RECIPE_DICT = dict(zip(list(map(lambda x: x.name, recipes)), recipes))

red_key = 'DRS_DATA_REDUC'


# =============================================================================
# Pre-load data for operations that take too much time
# =============================================================================
params = constants.load(__INSTRUMENT__)
red_index = tu.load_index_df(params[red_key])


# =============================================================================
# Define tests
# =============================================================================
tests = []

# -----------------------------------------------------------------------------
# Preprocessing Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Dark Master Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Badpixel Test
# -----------------------------------------------------------------------------
cal_badpix = DrsTest(drs_recipe=["cal_badpix_spirou.py"])
cal_badpix.name = "Bad Pixel Correction Recipe Test #1"
# cal_badpix.run_test = BadPixTest()
# TODO: Add other attributes here

tests.append(cal_badpix)

# -----------------------------------------------------------------------------
# Localisation Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Shape Master Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Shape (per night) Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Flat/Blaze Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Thermal Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Leak Master Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Wavelength Master Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Wavelgenth (per night) Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Leak (per night) Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Extraction Tests
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Make Telluric Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Fit Telluric Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Make Template Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CCF Test
# -----------------------------------------------------------------------------
