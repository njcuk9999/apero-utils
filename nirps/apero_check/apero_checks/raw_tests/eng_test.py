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
# Define classes
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)
class HdrKey:
    def __init__(self, name: str, header: Optional[str] = None,
                 dtype: Optional[str] = None, test: bool = True):
        self.name = str(name)

        if header is None:
            self.header_key = str(self.name)
        else:
            self.header_key = str(header)

        self.dtype = dtype
        self.in_header = False
        self.value = None
        self.test = test

    def read_header(self, hdr: fits.Header, filename: str):

        if self.header_key in hdr:
            raw_value = hdr[self.header_key]
            self.in_header = True
        else:
            raw_value = None
            self.in_header = False
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
            else:
                if raw_value is None:
                    self.value = ''
                else:
                    self.value = str(raw_value)
        except Exception as e:
            emsg = 'Error reading header key {0} with dtype {1} for file {2}'
            eargs = [self.header_key, self.dtype, filename]
            print(emsg.format(*eargs))
            print('Error was: {0}'.format(e))

    @staticmethod
    def print_keys(hdrlist: List['HdrKey']):
        # loop around header keys in provided list
        for hdrkey in hdrlist:
            # only show keys that we test
            if hdrkey.test:
                print(f'{hdrkey.name:30} {hdrkey.header_key}')


class EngTest:
    def __init__(self, name: str):
        self.name = name
        self.data = dict()
        self.vectors = dict()
        self.variables = dict()
        self.func = None
        self.fmsg = ''
        self.pmsg = ''
        self.calc = dict()
        self.active = True
        self.filters = dict()

    def run_test(self, tbl_dict: Dict[str, Any],
                 mask_dict: Dict[str, np.ndarray],
                 logger: List[str], passer: List[bool]
                 ) -> Tuple[List[str], List[bool]]:
        # if not active do not run test
        if not self.active:
            return logger, passer
        # log progress
        misc.log_msg(f'\t- Running sub-test: {self.name}', level='')
        # ---------------------------------------------------------------------
        # Step 1: Populate self.variables from tbl_dict
        # ---------------------------------------------------------------------
        # global mask (files that have both)
        global_mask = np.ones(len(tbl_dict['FILENAME'])).astype(bool)
        # Flag for no header entries
        not_found = False
        # loop around data and populate it from tbl_dict
        for vname in self.data:
            # get the key
            key = self.data[vname]
            # ------------------------------------------------------------------
            # some variables may not be in tbl_dict (just put these straight
            #  into self.variables
            if not isinstance(key, str):
                self.variables[vname] = self.data[vname]
                continue
            if key not in tbl_dict:
                self.variables[vname] = self.data[vname]
                continue
            # ------------------------------------------------------------------
            # get the value from the mask
            mask_value = mask_dict[key]
            # ------------------------------------------------------------------
            filter_str = ''
            # deal with filters
            for filter_name in self.filters:
                # get filter
                _filter = self.filters[filter_name]
                # if this filter is not in the table don't add it
                if filter_name not in tbl_dict:
                    continue
                # calculate filter mask
                filter_mask = np.in1d(tbl_dict[filter_name], _filter)
                # add to the string
                filter_str += f'{filter_name}=[{",".join(_filter)}] '
                # add this filter to the mask value
                mask_value &= filter_mask
            # ------------------------------------------------------------------
            # deal with no valid values
            if np.sum(mask_value) == 0:
                lmsg = f'No header entries for {key} (subtest={self.name}'
                lmsg += filter_str
                logger.append(lmsg)
                passer.append(True)
                not_found = True
            else:
                self.vectors[key] = np.array(tbl_dict[key][mask_value])
                self.variables[vname] = np.array(tbl_dict[key][mask_value])
            # push into global_mask
            global_mask &= mask_value
        # get a list of files that pass
        files = tbl_dict['FILENAME'][global_mask]
        # if we have no valid values for one of the parameters we stop here
        if not_found:
            return logger, passer
        # ---------------------------------------------------------------------
        # Step 2: Add calc variables to self.variables
        # ---------------------------------------------------------------------
        for vname in self.calc:
            # calculate the value
            value = self.calc[vname](**self.variables)
            # update self.variables
            self.variables[vname] = value

        # ---------------------------------------------------------------------
        # Step 3: Calculate the logic
        # ---------------------------------------------------------------------
        logic = self.func(**self.variables)
        # evaluate type (scalar logic should be boolean)
        if isinstance(logic, (bool, np.bool_)):
            vector = False
        else:
            vector = True
        # ---------------------------------------------------------------------
        # Step 4: Deal with pass/fail
        # ---------------------------------------------------------------------
        # deal with a scalar value
        if not vector:
            if logic:
                # append to logger
                logger.append(self.pmsg.format(**self.variables))
                # add pass to passer
                passer.append(True)
            else:
                # append to logger
                logger.append(self.fmsg.format(**self.variables))
                # add failure to passer
                passer.append(False)
        # deal with a vector value
        else:
            # pass if all are True
            if np.sum(logic) == len(logic):
                # append to logger
                logger.append(self.pmsg.format(**self.variables))
                # add pass to passer
                passer.append(True)
            else:
                # Add error to logger
                fmsg = self.fmsg.format(**self.variables)
                # add bad files to filelist
                for it, filename in enumerate(files):
                    if logic[it]:
                        continue
                    # add any vectors to the print out
                    vstr = ''
                    for vname in self.vectors:
                        vstr += f'\t{vname}={self.vectors[vname][it]}'
                    # push into the fail message
                    fmsg += f'\n\t {it+1}: {filename}' + vstr
                # append to logger
                logger.append(fmsg)
                # add failure to passer
                passer.append(False)
        # return logger and passer
        return logger, passer


# =============================================================================
# Define variables
# =============================================================================
# Set up header keys
HDR_KEYS = []
HDR_KEYS.append(HdrKey('DATE-OBS', dtype='time.iso', test=False))
HDR_KEYS.append(HdrKey('MJD-OBS', dtype='time.mjd', test=False))
HDR_KEYS.append(HdrKey('OBJECT', dtype='str', test=False))
HDR_KEYS.append(HdrKey('DPRTYPE',
                         header='HIERARCH ESO DPR TYPE',
                         dtype='str', test=False))
HDR_KEYS.append(HdrKey('TurboPumpStatus',
                         header='HIERARCH ESO INS SENS102 STAT',
                         dtype='bool'))
HDR_KEYS.append(HdrKey('EncloserTemperature',
                         header='HIERARCH ESO INS TEMP185 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('EncloserTemperatureSetpoint',
                         header='HIERARCH ESO INS TEMP187 VAL', dtype='float'))
HDR_KEYS.append(HdrKey('VacuumGauge1',
                         header='HIERARCH ESO INS PRES104 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('IsolationValve',
                         header='HIERARCH ESO INS SENS100 STAT',
                         dtype='bool'))
HDR_KEYS.append(HdrKey('Cryo1Status',
                         header='HIERARCH ESO INS SENS126',
                         dtype='bool'))
HDR_KEYS.append(HdrKey('Cryo2Status',
                         header='HIERARCH ESO INS SENS127',
                         dtype='bool'))
HDR_KEYS.append(HdrKey('WarningCryo1',
                         header='HIERARCH ESO INS SENS144 STAT',
                         dtype='str'))
HDR_KEYS.append(HdrKey('WarningCryo2',
                         header='HIERARCH ESO INS SENS146 STAT',
                         dtype='str'))
HDR_KEYS.append(HdrKey('FPtemperature_interior',
                         header='HIERARCH ESO INS TEMP14 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('FPtemperature_exterior',
                         header='HIERARCH ESO INS TEMP13 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('FPtemperature_setpoint',
                         header='HIERARCH ESO INS TEMP188 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('EncloserHeaterPower',
                         header='HIERARCH ESO INS SENS121 VAL',
                         dtype='float'))
HDR_KEYS.append(HdrKey('ScramblingStatus',
                       header='HIERARCH ESO INS2 AOS SCRAMB ST',
                       dtype='bool'))
HDR_KEYS.append(HdrKey('StretcherStatus',
                       header='HIERARCH ESO INS OPTI10 STAT',
                       dtype='str'))


# Set up engineering tests
# noinspection PyListCreation
ETESTS: Dict[str, EngTest] = dict()
# -----------------------------------------------------------------------------
ETESTS['etemp1'] = EngTest('test_enclosure_temperature')
ETESTS['etemp1'].data = dict(x='EncloserTemperature',
                             y='EncloserTemperatureSetpoint',
                             limit=0.1)
ETESTS['etemp1'].calc = dict(rms=lambda **k: np.nanstd(k['x']-k['y']))
ETESTS['etemp1'].func = lambda **k: k['rms'] < k['limit']
ETESTS['etemp1'].pmsg = 'RMS of enclosure temperature {rms:.1E} K < {limit:.1E} K'
ETESTS['etemp1'].fmsg = 'RMS of enclosure temperature {rms:.1E} K >= {limit:.1E} K'
# -----------------------------------------------------------------------------
ETESTS['etemp2'] = EngTest('test_enclosure_temperature_setpoint')
ETESTS['etemp2'].data = dict(x='EncloserTemperature',
                             y='EncloserTemperatureSetpoint',
                             limit=0.1)
ETESTS['etemp2'].calc = dict(mean_diff=lambda **k: np.nanmean(k['x'] - k['y']))
ETESTS['etemp2'].func = lambda **k: np.abs(k['mean_diff'] < k['limit'])
ETESTS['etemp2'].pmsg = 'Mean enclosure temperature diff {mean_diff:.1F} K < {limit:.1F} K'
ETESTS['etemp2'].fmsg = 'Mean enclosure temperature diff {mean_diff:.1F} K >= {limit:.1F} K'
# -----------------------------------------------------------------------------
ETESTS['vgaug1'] = EngTest('test_vacuum_gauge1')
ETESTS['vgaug1'].data = dict(x='VacuumGauge1',
                             limit=1.0e-4)
ETESTS['vgaug1'].calc = dict(maxx=lambda **k: np.nanmax(k['x']))
ETESTS['vgaug1'].func = lambda **k: k['x'] < k['limit']
ETESTS['vgaug1'].pmsg = 'max VacuumGaug1 {maxx:.2E} mbar < {limit:.2E} mbar'
ETESTS['vgaug1'].fmsg = 'max VacuumGaug1 {maxx:.2E} mbar >= {limit:.2E} mbar'
# -----------------------------------------------------------------------------
ETESTS['ivalv1'] = EngTest('test_isolation_valve')
ETESTS['ivalv1'].data = dict(x='IsolationValve')
ETESTS['ivalv1'].func = lambda **k: k['x'] == 0
ETESTS['ivalv1'].pmsg = 'Isolation valve closed'
ETESTS['ivalv1'].fmsg = 'Isolation valve open'
# -----------------------------------------------------------------------------
ETESTS['cryo1s'] = EngTest('test_cryo1_stat')
ETESTS['cryo1s'].data = dict(x='Cryo1Status')
ETESTS['cryo1s'].func = lambda **k: k['x'] == 0
ETESTS['cryo1s'].pmsg = 'Cryo1 on'
ETESTS['cryo1s'].fmsg = 'Cryo1 off'
# -----------------------------------------------------------------------------
ETESTS['cryo2s'] = EngTest('test_cryo2_stat')
ETESTS['cryo2s'].data = dict(x='Cryo2Status')
ETESTS['cryo2s'].func = lambda **k: k['x'] == 0
ETESTS['cryo2s'].pmsg = 'Cryo2 on'
ETESTS['cryo2s'].fmsg = 'Cryo2 off'
# -----------------------------------------------------------------------------
ETESTS['tpump1'] = EngTest('test_turbo_pump_status')
ETESTS['tpump1'].data = dict(x='TurboPumpStatus')
ETESTS['tpump1'].calc = dict(tsum=lambda **k: np.nansum(k['x']),
                             tmean=lambda **k: np.nanmean(k['x']),
                             tper=lambda **k: k['tmean'] * 100.0)
ETESTS['tpump1'].func = lambda **k: k['x'] == 0
ETESTS['tpump1'].pmsg = 'TurboPumpStatus always closed'
ETESTS['tpump1'].fmsg = 'TurboPumpStatus sometimes open ({tsum} times, {tper:.1f}% of time)'
# -----------------------------------------------------------------------------
ETESTS['cryo1w'] = EngTest('test_warning_cryo1')
ETESTS['cryo1w'].data = dict(x='WarningCryo1')
ETESTS['cryo1w'].calc = dict(tsum=lambda **k: np.nansum(k['x']))
ETESTS['cryo1w'].func = lambda **k: k['x'] == 0
ETESTS['cryo1w'].pmsg = 'WarningCryo1 always OK'
ETESTS['cryo1w'].fmsg = 'WarningCryo1 not empty ({0} times)'
ETESTS['cryo1w'].active = False
# -----------------------------------------------------------------------------
ETESTS['cryo2w'] = EngTest('test_warning_cryo2')
ETESTS['cryo2w'].data = dict(x='WarningCryo2')
ETESTS['cryo2w'].calc = dict(tsum=lambda **k: np.nansum(k['x']))
ETESTS['cryo2w'].func = lambda **k: k['x'] == 0
ETESTS['cryo2w'].pmsg = 'WarningCryo2 always OK'
ETESTS['cryo2w'].fmsg = 'WarningCryo2 not empty ({0} times)'
ETESTS['cryo2w'].active = False
# -----------------------------------------------------------------------------
ETESTS['fptemp'] = EngTest('test_fp_temperature')
ETESTS['fptemp'].data = dict(x='FPtemperature_interior',
                             limit=1.0e-2)
ETESTS['fptemp'].calc = dict(rms=lambda **k: np.nanstd(k['x']))
ETESTS['fptemp'].func = lambda **k: k['rms'] < k['limit']
ETESTS['fptemp'].pmsg = 'RMS FP temperature interior {rms:.2E} K < {limit:.2E} K)'
ETESTS['fptemp'].fmsg = 'RMS FP temperature interior {rms:.2E} K >= {limit:.2E} K)'
# -----------------------------------------------------------------------------
ETESTS['fptset'] = EngTest('test_fp_temperature_setpoint')
ETESTS['fptset'].data = dict(x='FPtemperature_interior',
                             y='FPtemperature_setpoint',
                             limit=0.005)
ETESTS['fptset'].calc = dict(rms=lambda **k: np.nanstd(k['x'] - k['y']))
ETESTS['fptset'].func = lambda **k: k['rms'] < k['limit']
ETESTS['fptset'].pmsg = 'RMS FP to setpoint {rms:.1E} K < {limit:.1E} K'
ETESTS['fptset'].fmsg = 'RMS FP to setpoint {rms:.1E} K >= {limit:.1E} K'
# -----------------------------------------------------------------------------
ETESTS['enhpow'] = EngTest('test_encloser_heater_power')
ETESTS['enhpow'].data = dict(x='EncloserHeaterPower',
                             limit=80.0)
ETESTS['enhpow'].calc = dict(maxx=lambda **k: np.nanmax(k['x']))
ETESTS['enhpow'].func = lambda **k: k['x'] < k['limit']
ETESTS['enhpow'].pmsg = 'Encloser heater power {maxx:.1f}% < {limit:.1f}%'
ETESTS['enhpow'].fmsg = 'Encloser heater power {maxx:.1f}% >= {limit:.1f}%'
# -----------------------------------------------------------------------------
ETESTS['scrsta'] = EngTest('ScramblingStatus')
ETESTS['scrsta'].data = dict(x='ScramblingStatus',
                             limit=True)
ETESTS['scrsta'].func = lambda **k: k['x'] == True
ETESTS['scrsta'].filters = dict(DPRTYPE=['OBJECT,FP', 'OBJECT,SKY', 'TELLURIC,SKY'])
ETESTS['scrsta'].pmsg = 'Scrambling Status Okay (=T)'
ETESTS['scrsta'].fmsg = 'Scrambling Status Not Okay (=F)'
# -----------------------------------------------------------------------------
ETESTS['scrdst'] = EngTest('StretcherStatus')
ETESTS['scrdst'].data = dict(x='StretcherStatus',
                             limit='ON')
ETESTS['scrdst'].calc = dict(cx=lambda **k: np.char.array(k['x']).strip())
ETESTS['scrdst'].func = lambda **k: k['cx'] == k['limit']
ETESTS['scrdst'].pmsg = 'Stretcher On'
ETESTS['scrdst'].fmsg = 'Stretcher Off'


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
    tbl_dict = dict(FILENAME=[])
    mask_dict = dict(FILENAME=[])
    # fill table looping through files
    for ifile in tqdm(range(len(files)), leave=False):
        # get the file
        tbl_dict['FILENAME'].append(files[ifile])
        mask_dict['FILENAME'].append(True)
        # get the header
        hdr = fits.getheader(files[ifile])
        # loop around keys and fill the table dictionary
        for it, key in enumerate(HDR_KEYS):
            key.read_header(hdr, files[ifile])
            # deal with new dictionary entry
            if key.name not in tbl_dict:
                tbl_dict[key.name] = []
                mask_dict[key.name] = []
            # push in value
            tbl_dict[key.name].append(key.value)
            mask_dict[key.name].append(key.in_header)
    # convert lists to numpy arrays
    for it, key in enumerate(tbl_dict):
        tbl_dict[key] = np.array(tbl_dict[key])
        mask_dict[key] = np.array(mask_dict[key]).astype(bool)
    # -------------------------------------------------------------------------
    # We add prints on the status of raw keywords:
    passer = []
    logger = []

    # -------------------------------------------------------------------------
    # Run tests
    # -------------------------------------------------------------------------
    # Run the engineering sub-tests
    for sub_test in ETESTS:
        # set up arguments for the sub-test functions
        test_args = [tbl_dict, mask_dict, logger, passer]
        # run test
        logger, passer = ETESTS[sub_test].run_test(*test_args)

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
