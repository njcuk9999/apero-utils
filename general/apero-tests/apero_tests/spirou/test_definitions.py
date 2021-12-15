"""
Definition of SPIRou recipe tests following APERO framework

@author: vandalt
"""
import apero_tests.global_subtests as gt
import apero_tests.utils as ut
from apero.core import constants
from apero_tests.drs_test import CACHEDIR, RECIPE_DICT, DrsTest

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

pp_key = "DRS_DATA_WORKING"
red_key = "DRS_DATA_REDUC"

# =============================================================================
# Pre-load data for operations that are relatively expensive
# =============================================================================
# ???: Maybe in future  we could have a single index df with extra index level
#  for parent dir (reduced, tmp, etc.)
params = constants.load(__INSTRUMENT__)
pp_log = ut.load_log_df(params[pp_key])
pp_index = ut.load_index_df(params[pp_key])
red_log = ut.load_log_df(params[red_key])
red_index = ut.load_index_df(params[red_key])

# If some files are on disk but not in index, we load their header and add them
# NOTE: For now this includes DEBUG files, which we all discard right after.
# We load them for consistency (because some DEBUG are already in index),
# but could add an option to skip them. Might matter for actual big datasets.
# Might matter less once we have SQL in 0.7 because can use query instead of
# file i/o

pp_missing_index, pp_not_found = ut.missing_index_headers(
    pp_index,
    instrument=__INSTRUMENT__,
    cache_dir=CACHEDIR,
    cache_suffix="_pp",
)
pp_full_index = ut.make_full_index(pp_index, pp_missing_index)

red_missing_index, red_not_found = ut.missing_index_headers(
    red_index,
    instrument=__INSTRUMENT__,
    cache_dir=CACHEDIR,
    cache_suffix="_red",
)
red_full_index = ut.make_full_index(red_index, red_missing_index)
master_calib_db = ut.load_db_entries(
    "CALIB", instrument=__INSTRUMENT__
)  # All calibdb entries
tellu_df = ut.load_db_entries(
    "TELLU", instrument=__INSTRUMENT__
)  # All telludb entries
# master_cdb_index = ut.load_db

# FUTURE: This should change in 0.7 with new way to match files per recipe
# Some global tests are done here (before per-recipe)
# - Report files not in index
# - Compare log and index with small report
# - Return index files that have a PID match in the logs
#   (the ones that can be processed at recipe level)
pp_index, pp_index_no_pid = ut.global_index_check(pp_full_index, pp_log)
red_index, red_index_no_pid = ut.global_index_check(red_full_index, red_log)
debug_mask = red_index.FILENAME.str.startswith("DEBUG")
red_index = red_index[~debug_mask]  # DEBUGs are not used in tests

full_index_no_pid = pp_index_no_pid.append(red_index_no_pid)
full_missing_index = pp_missing_index.append(red_missing_index)
full_not_found = pp_not_found.append(red_not_found)

# Get calib DB files used for files in index, and time difference
red_cdb_used_df = ut.get_cdb_df(
    red_index, params, cache_dir=CACHEDIR
)  # CDB keys and time difference with output files

# =============================================================================
# Define tests
# =============================================================================
tests = []

# -----------------------------------------------------------------------------
# "Global" test
# -----------------------------------------------------------------------------
# This one does not look like the others: it has no recipe and only has a few
# subtests
global_test = DrsTest(
    instrument=__INSTRUMENT__,
)
global_test.name = "Global test for all files"
global_test.test_id = "global_test"
global_test.subtest_list = [
    gt.GlobalIndexCheck(full_index_no_pid, global_test),
    gt.CheckMissingAdded(full_missing_index, global_test),
    gt.CheckNotFound(full_not_found, global_test),
]

tests.append(global_test)


# -----------------------------------------------------------------------------
# Preprocessing Test
# -----------------------------------------------------------------------------
cal_preprocess = DrsTest(
    drs_recipe=RECIPE_DICT["cal_preprocess_spirou.py"],
    pp_flag=True,
    all_log_df=pp_log,
    all_index_df=pp_index,
)

tests.append(cal_preprocess)


# -----------------------------------------------------------------------------
# Dark Master Test
# -----------------------------------------------------------------------------
cal_dark_master = DrsTest(
    drs_recipe=RECIPE_DICT["cal_dark_master_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_dark_master)


# -----------------------------------------------------------------------------
# Badpixel Test
# -----------------------------------------------------------------------------
cal_badpix = DrsTest(
    drs_recipe=RECIPE_DICT["cal_badpix_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)
# NOTE: Just an example subtest for future development (checks nothing)
cal_badpix.subtest_list.append(gt.CustomBadpixTest("This is a message "))

tests.append(cal_badpix)

# -----------------------------------------------------------------------------
# Localisation Test
# -----------------------------------------------------------------------------
cal_loc = DrsTest(
    drs_recipe=RECIPE_DICT["cal_loc_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
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
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_shape_master)

# -----------------------------------------------------------------------------
# Shape (per night) Test
# -----------------------------------------------------------------------------
cal_shape = DrsTest(
    drs_recipe=RECIPE_DICT["cal_shape_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_shape)

# -----------------------------------------------------------------------------
# Flat/Blaze Test
# -----------------------------------------------------------------------------
cal_flat = DrsTest(
    drs_recipe=RECIPE_DICT["cal_flat_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_flat)


# -----------------------------------------------------------------------------
# Thermal Test
# -----------------------------------------------------------------------------
cal_thermal = DrsTest(
    drs_recipe=RECIPE_DICT["cal_thermal_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
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
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_leak_master)


# -----------------------------------------------------------------------------
# Wavelength Master Test
# -----------------------------------------------------------------------------
cal_wave_master = DrsTest(
    drs_recipe=RECIPE_DICT["cal_wave_master_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_wave_master)


# -----------------------------------------------------------------------------
# Wavelgenth (per night) Test
# -----------------------------------------------------------------------------
cal_wave_night = DrsTest(
    drs_recipe=RECIPE_DICT["cal_wave_night_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_wave_night)


# -----------------------------------------------------------------------------
# Extraction Tests
# -----------------------------------------------------------------------------
cal_extract = DrsTest(
    drs_recipe=RECIPE_DICT["cal_extract_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_extract)


# -----------------------------------------------------------------------------
# Leak (per night) Test
# -----------------------------------------------------------------------------
cal_leak = DrsTest(
    drs_recipe=RECIPE_DICT["cal_leak_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

tests.append(cal_leak)


# -----------------------------------------------------------------------------
# Make Telluric Test
# -----------------------------------------------------------------------------

obj_mk_tellu = DrsTest(
    drs_recipe=RECIPE_DICT["obj_mk_tellu_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

# TODO: which OBJECTs was mk_tellu run on?
# was the correct template used (if used) should be that objects template
# how many have a template, how many don't?

tests.append(obj_mk_tellu)


# -----------------------------------------------------------------------------
# Fit Telluric Test
# -----------------------------------------------------------------------------

obj_fit_tellu = DrsTest(
    drs_recipe=RECIPE_DICT["obj_fit_tellu_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

# TODO: which OBJECTs was mk_tellu run on?
# was the correct template used (if used) should be that objects template
# how many have a template, how many don't?

tests.append(obj_fit_tellu)


# -----------------------------------------------------------------------------
# Make Template Test
# -----------------------------------------------------------------------------

obj_mk_template = DrsTest(
    drs_recipe=RECIPE_DICT["obj_mk_template_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

# TODO: how many objects were used in each template? does this match up to the
# number of tcorr files (output of fit_tellu)?

tests.append(obj_mk_template)


# -----------------------------------------------------------------------------
# CCF Test
# -----------------------------------------------------------------------------

cal_ccf = DrsTest(
    drs_recipe=RECIPE_DICT["cal_ccf_spirou.py"],
    all_log_df=red_log,
    all_index_df=red_index,
    all_master_calib_df=master_calib_db,
    all_tellu_df=tellu_df,
    all_cdb_used_df=red_cdb_used_df,
)

# TODO: what ccf step/ccf width/ccf target rv were used for each

tests.append(cal_ccf)
