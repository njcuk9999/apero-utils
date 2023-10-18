#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import copy
import glob
import os
from typing import Any, Dict

import numpy as np
from astropy.io import fits

from apero_checks.core import misc


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
    Test with observation directory exists

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # define parameters we use here
    raw_directory = params['raw dir']
    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    obsdir_path = os.path.join(raw_directory, obsdir)

    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {0} does not exist - TEST FAILED')
        return False

    files = glob.glob(os.path.join(obsdir_path, '*.fits'))

    if len(files) == 0:
        if log:
            # pass a True if no file is found on that night
            print('No files found for night {}'.format(obsdir))
        return False

    # create table to store keywords
    tbl = dict()

    # keywords to extract and table name
    keys = [['DATE-OBS', 'DATE-OBS'],
            ['MJD-OBS', 'MJD-OBS'],
            ['OBJECT', 'OBJECT'],
            ['HIERARCH ESO DPR TYPE', 'DPRTYPE'],
            ['ESO INS SENS102 STAT', 'TurboPumpStatus'],  # False
            ['ESO INS TEMP185 VAL', 'EncloserTemperature'],  # ~20 C
            ['ESO INS TEMP187 VAL', 'EncloserTemperatureSetpoint'],  # ~20 C
            ['ESO INS PRES104 VAL', 'VacuumGauge1'],  # < 1e-5 mbar
            ['ESO INS SENS100 STAT', 'IsolationValve'],  # False
            ['ESO INS SENS144 STAT', 'WarningCryo1'],  # ''
            ['ESO INS SENS146 STAT', 'WarningCryo2'],
            ['ESO INS TEMP14 VAL', 'FPtemperature_interior'],
            ['ESO INS TEMP13 VAL', 'FPtemperature_exterior'],
            ['ESO INS TEMP188 VAL', 'FPtemperature_setpoint']  # ''
            ]
    keys = np.array(keys)

    # define table columns
    h0 = fits.getheader(files[0])

    for i in range(len(keys[:, 0])):
        # loop on keys and which column to store them
        key = keys[i, 0]
        key2 = keys[i, 1]

        # if not present, we store a dummy value not to raise an error
        if key not in h0.keys():
            val = 0.0
        else:
            val = h0[key]

        # pad strings to avoid error if there are long strings in later files
        if type(val) == str:
            val = ' ' * 100
        tbl[key2] = np.array([val] * len(files))

    # fill table looping through files
    for ifile in range(len(files)):
        h = fits.getheader(files[ifile])

        for i in range(len(keys[:, 0])):
            key = keys[i, 0]  # keyword in header
            key2 = keys[i, 1]  # keyword in table
            # only if key is present in header
            if key in h.keys():
                tbl[key2][ifile] = copy.deepcopy(h[key])

    # change to True/False for strings
    for key in tbl.keys():
        if (tbl[key][0] == 'FALSE') or (tbl[key][0] == 'TRUE'):
            tbl[key] = tbl[key] == 'TRUE'

    # We add prints on the status of raw keywords:
    passed_log = []
    failed_log = []

    # on Enclosure Temperature
    ####################################################################################
    rms = np.nanstd(tbl['EncloserTemperature'] - tbl['EncloserTemperatureSetpoint'])
    # Enclosure temperature should be stable to 0.1 K rms
    qc_rms = 0.1
    qc = rms < qc_rms
    if qc:
        passed_log.append('\tRMS of enclosure temperature {:.1E} K  (<{:.1E} K)'.format(rms, qc_rms))
    else:
        failed_log.append('\tRMS of enclosure temperature {:.1E} K  (<{:.1E} K)'.format(rms, qc_rms))

    ####################################################################################
    # Enclosure temperature should be within 0.1 K of setpoint
    mean_diff = np.nanmean(tbl['EncloserTemperature'] - tbl['EncloserTemperatureSetpoint'])
    qc_mean = 0.1
    qc = np.abs(mean_diff) < qc_mean

    if qc:
        passed_log.append('\tmean enclosure temperature diff {:.1E} K  (<{:.1E} K)'.format(mean_diff, qc_mean))
    else:
        failed_log.append('\tmean enclosure temperature diff {:.1E} K  (<{:.1E} K)'.format(mean_diff, qc_mean))

    ####################################################################################
    # The VacuumGauge1 should see a really good vacuum otherwise there is a leak
    vacuum_gauge1 = np.nanmax(tbl['VacuumGauge1'])
    qc_vacuum_gauge1 = 1e-4
    qc = vacuum_gauge1 < qc_vacuum_gauge1

    if qc:
        passed_log.append('\tmax VacuumGauge1 {:.2E} mbar  (<{:.2E} mbar)'.format(vacuum_gauge1, qc_vacuum_gauge1))
    else:
        failed_log.append('\tmaxVacuumGauge1 {:.2E} mbar  (<{:.2E} mbar)'.format(vacuum_gauge1, qc_vacuum_gauge1))

    ####################################################################################
    # The IsolationValve should be closed, that's super bad if it's open
    valve_status = np.sum(tbl['IsolationValve'])

    if valve_status == 0:
        passed_log.append('\tIsolation valve always closed')
    else:
        failed_log.append('\tIsolation valve sometimes open ({} times)'.format(valve_status))

    ####################################################################################
    # The TurboPumpStatus should be 'off' on the fast majority of occurences (ask Philippe for details)
    turbo_status = np.sum(tbl['TurboPumpStatus'])

    if turbo_status == 0:
        passed_log.append('\tTurboPumpStatus always closed')
    else:
        failed_log.append('\tTurboPumpStatus sometimes open ({} times, '
                          '{:.1f}% of time)'.format(turbo_status, np.mean(tbl['TurboPumpStatus']) * 100.))

    ####################################################################################
    # Status is '' empty string when everything is OK
    # warning_cryo1_status = np.sum(tbl['WarningCryo1'] != '')
    #
    # if warning_cryo1_status == 0:
    #     passed_log.append('\tWarningCryo1 always OK')
    # else:
    #     failed_log.append('\tWarningCryo1 not empty ({} times)'.format(warning_cryo1_status))

    ####################################################################################
    # Status is '' empty string when everything is OK

    # warning_cryo2_status = np.sum(tbl['WarningCryo2'] != '')
    #
    # if warning_cryo2_status == 0:
    #     passed_log.append('\tWarningCryo2 always OK')
    # else:
    #     failed_log.append('\tWarningCryo2 not empty ({} times)'.format(warning_cryo2_status))

    ####################################################################################
    # FP temperature should be stable to 0.01 K rms
    rms_fp_temperature = np.nanstd(tbl['FPtemperature_interior'])
    qc_f_ptemperature = 1e-2
    qc = rms_fp_temperature < qc_f_ptemperature
    if qc:
        passed_log.append('\tRMS FP temperature interior {:.2E} K  (<{:.2E} K)'.format(rms_fp_temperature,
                                                                                       qc_f_ptemperature))
    else:
        failed_log.append('\tRMS FP temperature interior {:.2E} K  (<{:.2E} K)'.format(rms_fp_temperature,
                                                                                       qc_f_ptemperature))
    ####################################################################################
    # FP temperature should be within 0.005 K of setpoint
    rms = np.nanstd(tbl['FPtemperature_exterior'] - tbl['FPtemperature_setpoint'])
    qc_rms = 0.005
    qc = rms < qc_rms
    if qc:
        passed_log.append('\tRMS FP to setpoint {:.1E} K  (<{:.1E} K)'.format(rms, qc_rms))
    else:
        failed_log.append('\tRMS FP to setpoint {:.1E} K  (<{:.1E} K)'.format(rms, qc_rms))

    # construct a string for printing output
    passed_log = '\n'.join(passed_log)
    if len(passed_log) == 0:
        passed_log = '\tNone'
    failed_log = '\n'.join(failed_log)

    if len(failed_log) == 0:
        failed_log = '\tNone'
        passed_logical = True
    else:
        passed_logical = False

    if log:
        print('\n')
        print('**********************************************************************')
        print('QC for night {}'.format(obsdir))
        print('Passed QC:')
        print(passed_log)
        print('Failed QC:')
        print(failed_log)
        print('**********************************************************************')
        print('\n')

    # -------------------------------------------------------------------------
    return passed_logical


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
