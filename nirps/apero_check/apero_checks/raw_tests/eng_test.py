#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import glob
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from astropy.io import fits
from tqdm import tqdm

from apero_checks.core import misc


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)
class EngTest:
    def __init__(self, name: str, header: Optional[str] = None,
                 dtype: Optional[str] = None):

        self.name = str(name)

        if header is None:
            self.header_key = str(self.name)
        else:
            self.header_key = str(header)

        self.dtype = dtype
        self.in_header = False
        self.value = None

    def read_header(self, hdr: fits.Header, filename: str):

        if self.header_key in hdr:
            raw_value = hdr[self.header_key]
        else:
            raw_value = None
        # try to interpret the header key with the given dtype
        try:

            if self.dtype == 'time.iso':
                if raw_value is None:
                    self.value = ''
                else:
                    self.value = str(raw_value)
            elif self.dtype == 'time.mjd':
                if raw_value is None:
                    self.value = np.nan
                else:
                    self.value = float(raw_value)
            elif self.dtype == 'str':
                if raw_value is None:
                    self.value = ''
                else:
                    return str(raw_value)
            elif self.dtype == 'bool':
                if raw_value is None:
                    self.value = False
                elif raw_value in ['True', 'TRUE', 'true', 1, True]:
                    self.value = True
                else:
                    self.value = False
            elif self.dtype == 'float':
                if raw_value is None:
                    self.value = np.nan
                else:
                    self.value = float(raw_value)
        except Exception as e:
            emsg = 'Error reading header key {0} with dtype {1} for file {2}'
            eargs = [self.header_key, self.dtype, filename]
            print(emsg.format(*eargs))
            print('Error was: {0}'.format(e))


# Set up engineering header keys
# noinspection PyListCreation
ENG_TESTS = []
ENG_TESTS.append(EngTest('DATE-OBS', dtype='time.iso'))
ENG_TESTS.append(EngTest('MJD-OBS', dtype='time.mjd'))
ENG_TESTS.append(EngTest('OBJECT', dtype='str'))
ENG_TESTS.append(EngTest('DPRTYPE', header='HIERARCH ESO DPR TYPE',
                         dtype='str'))
ENG_TESTS.append(EngTest('TurboPumpStatus', header='ESO INS SENS102 STAT',
                         dtype='bool'))
ENG_TESTS.append(EngTest('EncloserTemperature', header='ESO INS TEMP185 VAL',
                         dtype='float'))
ENG_TESTS.append(EngTest('EncloserTemperatureSetpoint',
                         header='ESO INS TEMP187 VAL', dtype='float'))
ENG_TESTS.append(EngTest('VacuumGauge1', header='ESO INS PRES104 VAL',
                         dtype='float'))
ENG_TESTS.append(EngTest('IsolationValve', header='ESO INS SENS100 STAT',
                         dtype='bool'))
ENG_TESTS.append(EngTest('WarningCryo1', header='ESO INS SENS144 STAT',
                         dtype='str'))
ENG_TESTS.append(EngTest('WarningCryo2', header='ESO INS SENS146 STAT',
                         dtype='str'))
ENG_TESTS.append(EngTest('FPtemperature_interior', header='ESO INS TEMP14 VAL',
                         dtype='float'))
ENG_TESTS.append(EngTest('FPtemperature_exterior', header='ESO INS TEMP13 VAL',
                         dtype='float'))
ENG_TESTS.append(EngTest('FPtemperature_setpoint', header='ESO INS TEMP188 VAL',
                         dtype='float'))
ENG_TESTS.append(EngTest('EncloserHeaterPower', header='ESO INS SENS121 VAL',
                         dtype='float'))


# =============================================================================
# Define test functions
# =============================================================================
def test_enclosure_temperature(tbl_dict: Dict[str, np.ndarray],
                               mask_dict: Dict[str, np.ndarray],
                               logger: List[str], passer: List[bool]
                               ) -> Tuple[List[str], List[bool]]:
    """
    Test Enclosure Temperature

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'EncloserTemperature'
    key2 = 'EncloserTemperatureSetpoint'

    rms = np.nanstd(tbl_dict[key1][mask_dict[key1]] -
                    tbl_dict[key2][mask_dict[key1]])
    # Enclosure temperature should be stable to 0.1 K rms
    qc_rms = 0.1
    qc_passed = rms < qc_rms

    msg = 'RMS of enclosure temperature {0:.1E} K  ({1}{2:.1E} K)'
    if qc_passed:
        margs = [rms, '<', qc_rms]

    else:
        margs = [rms, '>', qc_rms]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_enclosure_temperature_setpoint(tbl_dict: Dict[str, Any],
                                        mask_dict: Dict[str, np.ndarray],
                                        logger: List[str], passer: List[bool]
                                        ) -> Tuple[List[str], List[bool]]:
    """
    Test Enclosure temperature should be within 0.1 K of setpoint

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'EncloserTemperature'
    key2 = 'EncloserTemperatureSetpoint'

    mean_diff = np.nanmean(tbl_dict[key1][mask_dict[key1]] -
                           tbl_dict[key2][mask_dict[key1]])

    qc_mean = 0.1
    qc_passed = np.abs(mean_diff) < qc_mean

    msg = 'mean enclosure temperature diff {0:.1E} K  ({1}{2:.1E} K)'
    if qc_passed:
        margs = [mean_diff, '<', qc_mean]
    else:
        margs = [mean_diff, '>', qc_mean]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_vacuum_gauge1(tbl_dict: Dict[str, Any], logger: List[str],
                       mask_dict: Dict[str, np.ndarray],
                       passer: List[bool]) -> Tuple[List[str], List[bool]]:
    """
    Test the VacuumGauge1 should see a really good vacuum otherwise there
    is a leak

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'VacuumGauge1'

    vacuum_gauge1 = np.nanmax(tbl_dict[key1][mask_dict[key1]])

    qc_vacuum_gauge1 = 1e-4
    qc_passed = vacuum_gauge1 < qc_vacuum_gauge1

    msg = 'max VacuumGauge1 {0:.2E} mbar  ({1}{2:.2E} mbar)'
    if qc_passed:
        margs = [vacuum_gauge1, '<', qc_vacuum_gauge1]
    else:
        margs = [vacuum_gauge1, '>', qc_vacuum_gauge1]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_isolation_valve(tbl_dict: Dict[str, Any], logger: List[str],
                         mask_dict: Dict[str, np.ndarray],
                         passer: List[bool]) -> Tuple[List[str], List[bool]]:
    """
    The IsolationValve should be closed, that's super bad if it's open

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'IsolationValve'

    isolation_valve = np.nansum(tbl_dict[key1][mask_dict[key1]])

    qc_passed = isolation_valve == 0
    if qc_passed:
        msg = 'Isolation valve always closed'
        margs = []
    else:
        msg = 'Isolation valve sometimes open ({0} times)'
        margs = [isolation_valve]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_turbo_pump_status(tbl_dict: Dict[str, Any], logger: List[str],
                           mask_dict: Dict[str, np.ndarray],
                           passer: List[bool]) -> Tuple[List[str], List[bool]]:
    """
    The TurboPumpStatus should be 'off' on the fast majority of occurences
    (ask Philippe for details)

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'TurboPumpStatus'

    turbo_status = np.nansum(tbl_dict[key1][mask_dict[key1]])

    qc_passed = turbo_status == 0
    if qc_passed:
        msg = 'TurboPumpStatus always closed'
        margs = []
    else:
        msg = 'TurboPumpStatus sometimes open ({0} times, {1:.1f}% of time)'
        margs = [turbo_status, np.nanmean(tbl_dict['TurboPumpStatus']) * 100.]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_warning_cryo1(tbl_dict: Dict[str, Any], logger: List[str],
                       mask_dict: Dict[str, np.ndarray],
                       passer: List[bool]) -> Tuple[List[str], List[bool]]:
    """
    Test the Cryo1 warning

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'WarningCryo1'

    warning_cryo1_status = np.nansum(tbl_dict[key1][mask_dict[key1]] != '')

    qc_passed = warning_cryo1_status == 0
    if qc_passed:
        msg = 'WarningCryo1 always OK'
        margs = []
    else:
        msg = 'WarningCryo1 not empty ({0} times)'
        margs = [warning_cryo1_status]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_warning_cryo2(tbl_dict: Dict[str, Any],
                       mask_dict: Dict[str, np.ndarray],
                       logger: List[str], passer: List[bool]
                       ) -> Tuple[List[str], List[bool]]:
    """
    Test the Cryo2 warning

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'WarningCryo2'

    warning_cryo2_status = np.nansum(tbl_dict[key1][mask_dict[key1]] != '')

    qc_passed = warning_cryo2_status == 0
    if qc_passed:
        msg = 'WarningCryo2 always OK'
        margs = []
    else:
        msg = 'WarningCryo2 not empty ({0} times)'
        margs = [warning_cryo2_status]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_fp_temperature(tbl_dict: Dict[str, Any],
                        mask_dict: Dict[str, np.ndarray],
                        logger: List[str],
                        passer: List[bool]) -> Tuple[List[str], List[bool]]:
    """
    Test FP temperature should be stable to 0.01 K rms

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'FPtemperature_interior'

    rms_fp_temperature = np.nanstd(tbl_dict[key1][mask_dict[key1]])

    qc_f_ptemperature = 1e-2
    qc_passed = rms_fp_temperature < qc_f_ptemperature
    if qc_passed:
        msg = 'RMS FP temperature interior {0:.2E} K  ({1}{2:.2E} K)'
        margs = [rms_fp_temperature, '<', qc_f_ptemperature]
    else:
        msg = 'RMS FP temperature interior {0:.2E} K  ({1}{2:.2E} K)'
        margs = [rms_fp_temperature, '>', qc_f_ptemperature]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_fp_temperature_setpoint(tbl_dict: Dict[str, Any],
                                 mask_dict: Dict[str, np.ndarray],
                                 logger: List[str], passer: List[bool]
                                 ) -> Tuple[List[str], List[bool]]:
    """
    FP temperature should be within 0.005 K of setpoint

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'FPtemperature_interior'
    key2 = 'FPtemperature_setpoint'

    rms = np.nanstd(tbl_dict[key1][mask_dict[key1]] -
                    tbl_dict[key2][mask_dict[key1]])

    qc_rms = 0.005
    qc_passed = rms < qc_rms
    if qc_passed:
        msg = 'RMS FP to setpoint {0:.1E} K  ({1}{2:.1E} K)'
        margs = [rms, '<', qc_rms]
    else:
        msg = 'RMS FP to setpoint {0:.1E} K  ({1}{2:.1E} K)'
        margs = [rms, '>', qc_rms]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


def test_encloser_heater_power(tbl_dict: Dict[str, Any],
                               mask_dict: Dict[str, np.ndarray],
                               logger: List[str], passer: List[bool]
                               ) -> Tuple[List[str], List[bool]]:
    """
    Test Encloser heater power should be less than 80%

    :param tbl_dict: dictionary of table values from headers of all files
    :param logger: list of strings with reason for passing or failing
    :param passer: list of bools, whether a test passed (True) or failed (False)

    :return: Tuple, the updated 1. logger and 2. passer
    """
    key1 = 'EncloserHeaterPower'

    max_heater_power = np.nanmax(tbl_dict[key1][mask_dict[key1]])

    heater_power_thres = 80.0
    qc_passed = max_heater_power < heater_power_thres
    if qc_passed:
        msg = 'Encloser heater power {0:.1f}% ({1}{2:.1f}%)'
        margs = [max_heater_power, '<', heater_power_thres]
    else:
        msg = 'Encloser heater power {0:.1f}% ({1}{2:.1f}%)'
        margs = [max_heater_power, '>', heater_power_thres]
    # add to the logger
    logger.append(msg.format(*margs))
    # add to the passed
    passer.append(qc_passed)
    # return logger and passer
    return logger, passer


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
    # create the observation directory path
    obsdir_path = os.path.join(raw_directory, obsdir)
    # deal with no observation directory
    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {0} does not exist - skipping test')
        return True
    # get a list of files in this observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))
    # deal with no files
    if len(files) == 0:
        if log:
            # pass a True if no file is found on that night
            print('No files found for night {} - skipping test'.format(obsdir))
        return True
    # -------------------------------------------------------------------------
    # create table to store keywords
    tbl_dict = dict()
    mask_dict = dict()
    # fill table looping through files
    for ifile in tqdm(range(len(files))):
        # get the header
        hdr = fits.getheader(files[ifile])
        # loop around keys and fill the table dictionary
        for it, key in enumerate(ENG_TESTS):
            key.read_header(hdr, files[ifile])
            if key.in_header:
                if key.name not in tbl_dict:
                    tbl_dict[key.name] = []
                    mask_dict[key.name] = []
                tbl_dict[key.name].append(key.value)
                mask_dict[key.name].append(True)
            else:
                if key.name not in tbl_dict:
                    tbl_dict[key.name] = []
                    mask_dict[key.name] = []
                tbl_dict[key.name].append(np.nan)
                mask_dict[key.name].append(False)

    # convert lists to numpy arrays
    for it, key in enumerate(ENG_TESTS):
        tbl_dict[key] = np.array(tbl_dict[key])
        mask_dict[key] = np.array(mask_dict[key]).astype(bool)

    # -------------------------------------------------------------------------
    # We add prints on the status of raw keywords:
    passer = []
    logger = []

    # -------------------------------------------------------------------------
    # Test Enclosure Temperature
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_enclosure_temperature(*test_args)

    # -------------------------------------------------------------------------
    # Test Enclosure temperature should be within 0.1 K of setpoint
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_enclosure_temperature(*test_args)

    # -------------------------------------------------------------------------
    # Test the VacuumGauge1 should see a really good vacuum otherwise there
    # is a leak
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_vacuum_gauge1(*test_args)

    # -------------------------------------------------------------------------
    # The IsolationValve should be closed, that's super bad if it's open
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_isolation_valve(*test_args)

    # -------------------------------------------------------------------------
    # The TurboPumpStatus should be 'off' on the fast majority of occurences
    # (ask Philippe for details)
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_turbo_pump_status(*test_args)

    # -------------------------------------------------------------------------
    # Test the Cryo1 warning
    # -------------------------------------------------------------------------
    # # set up arguments for the sub-test functions
    # test_args = [tbl_dict, mask_dict, logger, passer]
    # # run test
    # logger, passer = test_warning_cryo1(*test_args)

    # -------------------------------------------------------------------------
    # Test the Cryo2 warning
    # -------------------------------------------------------------------------
    # # set up arguments for the sub-test functions
    # test_args = [tbl_dict, mask_dict, logger, passer]
    # # run test
    # logger, passer = test_warning_cryo2(*test_args)

    # -------------------------------------------------------------------------
    # Test FP temperature should be stable to 0.01 K rms
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_fp_temperature(*test_args)

    # -------------------------------------------------------------------------
    # FP temperature should be within 0.005 K of setpoint
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_fp_temperature_setpoint(*test_args)

    # -------------------------------------------------------------------------
    # Test Encloser heater power should be less than 80%
    # -------------------------------------------------------------------------
    # set up arguments for the sub-test functions
    test_args = [tbl_dict, mask_dict, logger, passer]
    # run test
    logger, passer = test_encloser_heater_power(*test_args)

    # -------------------------------------------------------------------------
    # construct a string for printing output
    # -------------------------------------------------------------------------
    # storage for combined log
    passed_log = []
    failed_log = []
    # loop around tests
    for it in range(len(passer)):
        if passer[it]:
            passed_log.append('\t' + logger[it])
        else:
            failed_log.append('\t' + logger[it])
    # deal with return
    if len(failed_log) == 0:
        passed_logical = True
    else:
        passed_logical = False
    # -------------------------------------------------------------------------
    if log:
        print('\n')
        print('*' * 50)
        print('QC for night {}'.format(obsdir))
        print('Passed QC:')
        print('\n'.join(passed_log))
        print('Failed QC:')
        print('\n'.join(failed_log))
        print('*' * 50)
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
