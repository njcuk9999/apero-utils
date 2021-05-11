"""
Definition of SPIRou recipe tests following APERO framework

@author: vandalt
"""
import apero_tests.utils as ut
from apero.core import constants
from apero.core.instruments.spirou.recipe_definitions import recipes
from apero_tests.drs_test import CACHEDIR, DrsTest

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
__NAME__ = "apero_tests.spirou.test_definitions.py"
__INSTRUMENT__ = "SPIROU"
# Get constants
Constants = constants.load(__INSTRUMENT__)
# Get version and author
__version__ = Constants["DRS_VERSION"]
__author__ = Constants["AUTHORS"]
__date__ = Constants["DRS_DATE"]
__release__ = Constants["DRS_RELEASE"]

RECIPE_DICT = dict(zip(list(map(lambda x: x.name, recipes)), recipes))

red_key = "DRS_DATA_REDUC"

# =============================================================================
# Pre-load data for operations that are relatively expensive
# =============================================================================
# ???: Maybe in future  we could have a single index df with extra index level
#  for parent dir (reduced, tmp, etc.)
# TODO: Better handle non-reduced log/index stuff (see ??? above)
params = constants.load(__INSTRUMENT__)
red_log = ut.load_log_df(params[red_key])
red_index = ut.load_index_df(params[red_key])

# If some files are on disk but not in index, we load their header and add them
# NOTE: For now this includes DEBUG files, which we all discard right after.
# We load them for consistency (because some DEBUG are already in index),
# but could add an option to skip them
red_missing_index = ut.missing_index_headers(
    red_index, instrument=__INSTRUMENT__, cache_dir=CACHEDIR
)
red_full_index = ut.make_full_index(red_index, red_missing_index)
master_calib_db = ut.load_db("CALIB", instrument=__INSTRUMENT__)  # All calibdb entries

# Some global tests are done here (before per-recipe)
# NOTE: This should change in 0.7. Will have new way to match files per recipe
# - Report files not in index
# - Compare log and index with small report
# - Return index files that have a PID match in the logs
#   (the ones that can be processed at index level)
red_index = ut.global_index_check(red_full_index, red_log)
debug_mask = red_index.FILENAME.str.startswith("DEBUG")
red_index = red_index[~debug_mask]  # DEBUGs are not used in testes
red_cdb_used_df = ut.get_cdb_df(
    red_index, params, cache_dir=CACHEDIR
)  # CDB keys and time difference with output files

# =============================================================================
# Define tests
# =============================================================================
tests = []

# TODO: Make compatible with pp and dark master testse
# -----------------------------------------------------------------------------
# Preprocessing Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Dark Master Test
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Badpixel Test
# -----------------------------------------------------------------------------
cal_badpix = DrsTest(
    drs_recipe=RECIPE_DICT["cal_badpix_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_badpix)

# -----------------------------------------------------------------------------
# Localisation Test
# -----------------------------------------------------------------------------
cal_loc = DrsTest(
    drs_recipe=RECIPE_DICT["cal_loc_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_loc)

# -----------------------------------------------------------------------------
# Shape Master Test
# -----------------------------------------------------------------------------
cal_shape_master = DrsTest(
    drs_recipe=RECIPE_DICT["cal_shape_master_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_shape_master)

# -----------------------------------------------------------------------------
# Shape (per night) Test
# -----------------------------------------------------------------------------
cal_shape = DrsTest(
    drs_recipe=RECIPE_DICT["cal_shape_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_shape)

# -----------------------------------------------------------------------------
# Flat/Blaze Test
# -----------------------------------------------------------------------------
cal_flat = DrsTest(
    drs_recipe=RECIPE_DICT["cal_flat_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_flat)


# -----------------------------------------------------------------------------
# Thermal Test
# -----------------------------------------------------------------------------
cal_thermal = DrsTest(
    drs_recipe=RECIPE_DICT["cal_thermal_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_thermal)


# -----------------------------------------------------------------------------
# Leak Master Test
# -----------------------------------------------------------------------------
cal_leak_master = DrsTest(
    drs_recipe=RECIPE_DICT["cal_leak_master_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_leak_master)


# -----------------------------------------------------------------------------
# Wavelength Master Test
# -----------------------------------------------------------------------------
cal_wave_master = DrsTest(
    drs_recipe=RECIPE_DICT["cal_wave_master_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_wave_master)


# -----------------------------------------------------------------------------
# Wavelgenth (per night) Test
# -----------------------------------------------------------------------------
cal_wave_night = DrsTest(
    drs_recipe=RECIPE_DICT["cal_wave_night_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_wave_night)


# -----------------------------------------------------------------------------
# Leak (per night) Test
# -----------------------------------------------------------------------------
cal_leak = DrsTest(
    drs_recipe=RECIPE_DICT["cal_leak_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_leak)


# -----------------------------------------------------------------------------
# Extraction Tests
# -----------------------------------------------------------------------------
cal_extract = DrsTest(
    drs_recipe=RECIPE_DICT["cal_extract_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_cdb_used_df=red_cdb_used_df,
)
# TODO: Add other attributes here (after structure is more complete)
# cal_badpix.run_test = BadPixTest()

tests.append(cal_extract)


# TODO: Add compat for science recipes
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
