#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict, List, Tuple

import numpy as np
import git

from apero_checks.core import base
from apero_checks.core import red_functions

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================
def calib_check(params: Any, recipe: Any, obsdir: str, log: bool = False
                ) -> Tuple[Dict[str, dict], Dict[str, list], List[str]]:
    """
    Code below is similar to apero.tools.module.processing.drs_precheck.
    file_check function

    :param params: the APERO parameter dictionary
    :param recipe: the APERO recipe
    :param obsdir: str, the observation directory
    :param log: bool, whether to print message to screen

    :return: tuple, 1. the calibration count for each obsdir
                    2. the calibration times for each obsdir
                    3. a list of bad calibration nights
    """
    # imports after apero profile update
    from apero.core.core import drs_database
    from apero.core.core import drs_text
    from apero.science import preprocessing as prep
    from apero.tools.module.processing import drs_processing
    from apero.tools.module.processing import drs_precheck
    from apero.science import telluric
    # -------------------------------------------------------------------------
    # construct the index database instance
    findexdbm = drs_database.FileIndexDatabase(params)
    findexdbm.load_db()
    # -------------------------------------------------------------------------
    # get odometer reject list (if required)
    # -------------------------------------------------------------------------
    # get whether the user wants to use reject list
    _use_odo_reject = params['USE_REJECTLIST']
    # get the odometer reject list
    odo_reject_list = []
    if not drs_text.null_text(_use_odo_reject, ['', 'None']):
        if drs_text.true_text(_use_odo_reject):
            odo_reject_list = prep.get_file_reject_list(params)
    # -------------------------------------------------------------------------
    # get the conditions based on params
    # -------------------------------------------------------------------------
    # TODO: Remove this try once online is on version v0.7.287+
    try:
        condition, _ = drs_processing.gen_global_condition(params, findexdbm,
                                                           odo_reject_list,
                                                           log=log)
    except Exception:
        condition = drs_processing.gen_global_condition(params, findexdbm,
                                                        odo_reject_list)
    # -------------------------------------------------------------------------
    # get telluric stars and non-telluric stars
    # -------------------------------------------------------------------------
    # TODO: Remove this try once online is on version v0.7.287+
    try:
        # get a list of all objects from the file index database
        all_objects = drs_processing.get_uobjs_from_findex(params, findexdbm)
        # get all telluric stars
        tstars = telluric.get_tellu_include_list(params, all_objects=all_objects)
        # get all other stars
        ostars = drs_processing.get_non_telluric_stars(params, all_objects, tstars)
    except Exception:
        # get all telluric stars
        tstars = telluric.get_tellu_include_list(params)
        # get all other stars
        ostars = drs_processing.get_non_telluric_stars(params, findexdbm, tstars)
    # -------------------------------------------------------------------------
    # get a list of calibration files
    # -------------------------------------------------------------------------
    # set up the obsdirs
    obsdirs = np.array([obsdir])
    # return the calib check fucntion from apero
    return drs_precheck.calib_check(params, recipe, tstars,  ostars, obsdirs,
                                    condition, findexdbm, log=log)


def test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Test with observation directory exists

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """

    # get the apero profiles run.ini file
    runfile = params['processing']['run file']
    # deal with no run file
    if runfile is None:
        emsg = ('APERO_CALIB_TEST error: "processing.run file" must be '
                'defined in yaml')
        raise base.AperoChecksError(emsg)
    # update apero profile
    apero_params = red_functions.update_apero_profile(params)
    # get the proxy apero recipe
    apero_recipe = red_functions.get_apero_proxy_recipe(apero_params)
    # update apero params with parameters that normally come from run.ini file
    apero_params = red_functions.add_run_ini_params(apero_params,
                                                    apero_recipe, runfile)
    # do the apero calib file check
    cout = calib_check(apero_params, apero_recipe, obsdir, log=log)
    calib_count, calib_times, bad_calib_nights = cout
    # -------------------------------------------------------------------------
    # we use bad_calib_nights to register pass/fail
    if len(bad_calib_nights) == 0:
        return True
    else:
        return False



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-07-01'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
