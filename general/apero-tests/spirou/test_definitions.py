"""
Definition of SPIRou recipe tests following APERO framework

@author: vandalt
"""
from apero.core import constants
from apero.core.instruments.spirou.recipe_definitions import recipes

# TODO: Make sure this works for dev
from ..drs_test import DrsTest
from .. import test_utils as ut

# Something similar?
# from .preprocessing import PPTest
# from .darkmaster import DarkMTest

# from .localisation import LocTest
# from .shapemaster import ShapeMTest
# from .shape import ShapeTest
# from .flat import FlatTest
# from .thermal import ThermalTest
# from .leakmaster import LeakMTest
# from .wavelengthmaster import WaveMTest
# from .wavelength import WaveTest

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
# Pre-load data for operations that are relatively expensive
# =============================================================================
# ???: Maybe in future  we could have a single index df with extra index level
#  for parent dir (reduced, tmp, etc.)
params = constants.load(__INSTRUMENT__)
red_log = ut.load_log_df(params[red_key])
red_index = ut.load_index_df(params[red_key])
red_missing_index = ut.missing_index_headers(red_index,
                                             instrument=__INSTRUMENT__)
red_full_index = ut.make_full_index(red_index, red_missing_index)
master_calib_db = ut.load_db("CALIB", instrument=__INSTRUMENT__)

# TODO: Before per-recipe tests, we have some cleanups to do, could either be a general
# DrsTest object or just few lines/util function here
# - Report files not in index
# - Comapre log and index with small report

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
# TODO: Add other attributes here (after structure is more complete)

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
