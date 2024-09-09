#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
from typing import Any, Dict

from apero_checks.core import base
from apero_checks.core import apero_functions

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
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
    # get the apero profiles run.ini file
    runfile = params['processing']['run file']
    # deal with no run file
    if runfile is None:
        emsg = ('APERO_CALIB_TEST error: "processing.run file" must be '
                'defined in yaml')
        raise base.AperoChecksError(emsg)
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
        # all print out messages must be wrapped in if log
        if log:
            # print error
            print('\nSome objects must be added to astrometric database '
                  '(via apero_astrometrics) or added to the reject list.\n')
            # loop around the rows
            for row in range(len(unfound_table)):
                # print the object
                msg = ('\t{0}\t{1:30s}\t(APERO: {2})'
                       '\tLAST[{3}, {4}, {5}]')
                margs = [row + 1,
                         ' or '.join(unfound_table['Original Names'][row]),
                         unfound_table['Apero Name'][row],
                         unfound_table['Last Run ID'][row],
                         unfound_table['Last PI Name'] [row],
                         unfound_table['Last Obs Date'][row]]
                # print the message
                print(msg.format(*margs))
        # return False
        return False
    else:
        # all print out messages must be wrapped in if log
        if log:
            print('No unfound objects.')
        # return True
        return True


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
