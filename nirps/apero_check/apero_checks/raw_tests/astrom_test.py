#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict, Optional, Tuple

from apero_checks.core import base
from apero_checks.core import apero_functions

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)

# Cache values, we only need to run test once (as it affects all nights)
ASTROM_CACHE: Optional[Tuple[bool, str]] = None

# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> Tuple[bool, str]:
    """
    Blank test - this tests whether test was run (should always return True)
    All other tests should return True or False, and only print messages if
    log is True.

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    global ASTROM_CACHE
    # get the apero profiles run.ini file
    runfile = params['processing']['run file']
    # deal with no run file
    if runfile is None:
        emsg = ('APERO_CALIB_TEST error: "processing.run file" must be '
                'defined in yaml')
        raise base.AperoChecksError(emsg)
    # deal with result cached from before (this can be done as this test does
    #  not depend on obsdir or any input from the user)
    if ASTROM_CACHE is not None:
        return ASTROM_CACHE[0], ASTROM_CACHE[1]

    # update apero profile
    apero_params = apero_functions.update_apero_profile(params)
    # get the proxy apero recipe
    apero_recipe = apero_functions.get_apero_proxy_recipe(apero_params)
    # update apero params with parameters that normally come from run.ini file
    apero_params = apero_functions.add_run_ini_params(apero_params,
                                                      apero_recipe, runfile)
    # we do not use obsdir here
    _ = obsdir
    # imports after apero profile update
    from apero.core.core import drs_database
    from apero.tools.module.processing import drs_precheck
    from apero.tools.module.database import manage_databases
    # -------------------------------------------------------------------------
    # update database
    manage_databases.update_object_database(apero_params, log=False)
    # -------------------------------------------------------------------------
    # construct the index database instance
    findexdbm = drs_database.FileIndexDatabase(apero_params)
    findexdbm.load_db()
    # get the unfound table
    unfound_table = drs_precheck.obj_check(apero_params, findexdbm, log=False)
    # -------------------------------------------------------------------------
    # We get a False condition if there are unfound objects
    if len(unfound_table) > 0:
        # print error
        out_msg = ('\nSome objects must be added to astrometric database '
              '(via apero_astrometrics) or added to the reject list.\n')
        # loop around the rows
        for row in range(len(unfound_table)):
            # print the object
            msg = ('\n\t{0}\t{1:30s}\t(APERO: {2})'
                   '\tLAST[{3}, {4}, {5}]')
            margs = [row + 1,
                     ' or '.join(unfound_table['Original Names'][row]),
                     unfound_table['Apero Name'][row],
                     unfound_table['Last Run ID'][row],
                     unfound_table['Last PI Name'] [row],
                     unfound_table['Last Obs Date'][row]]
            # print the message
            out_msg += (msg.format(*margs))
        # all print out messages must be wrapped in if log
        if log:
            print(out_msg)
        # cache results
        ASTROM_CACHE = [False, out_msg]
        # return False
        return False, out_msg
    else:
        # all print out messages must be wrapped in if log
        out_msg = 'No unfound objects.'
        if log:
            print(out_msg)
        # cache results
        ASTROM_CACHE = [True, out_msg]
        # return True
        return True, out_msg


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-03-15'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
