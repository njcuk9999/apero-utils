#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-02 at 11:47

@author: cook
"""
import os
import shutil
from typing import List

import numpy as np
from tqdm import tqdm

from apero.base import base
from apero.core import constants
from apero.core.core import drs_database
from apero.core.core import drs_log

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'minidata_selector.py'
__INSTRUMENT__ = 'None'
__PACKAGE__ = base.__PACKAGE__
__version__ = base.__version__
__author__ = base.__author__
__date__ = base.__date__
__release__ = base.__release__
# Get Logging function
WLOG = drs_log.wlog
ParamDict = constants.ParamDict
# whether we are in test mode
TEST = False
# Define the required science targets for a mini dataset
SCIENCE_TARGETS = dict()
SCIENCE_TARGETS['SPIROU'] = ['GL699']
SCIENCE_TARGETS['NIRPS_HE'] = ['PROXIMA']
SCIENCE_TARGETS['NIRPS_HA'] = ['PROXIMA']
# Directory to copy files to
OUTPUTDIR = dict()
OUTPUTDIR['SPIROU'] = '/scratch2/spirou/drs-data/common/minidata2_xxs'
OUTPUTDIR['NIRPS_HE'] = '/cosmos99/nirps/apero-data/misc/common/nirps_he_xxs'
OUTPUTDIR['NIRPS_HA'] = '/cosmos99/nirps/apero-data/misc/common/nirps_ha_xxs'
# -----------------------------------------------------------------------------
# define a list of hot stars
TELLURIC_TARGETS = dict()
TELLURIC_TARGETS['SPIROU'] = [
    "74PscB", "31Cas", "gamTri", "HR875", "HR1314", "pi.02Ori", "HR1832", "zetLep",
    "HR2180", "HR2209", "24Lyn", "HR3131", "33Lyn", "etaPyx", "23LMi", "lLeo", "phiLeo",
    "HR4687", "HR4722", "zetVir", "82UMa", "HD130917", "betSer", "HR6025", "HD159170",
    "gamSct", "51Dra", "iotCyg", "omiCapA", "chiCap", "17Peg", "HR8489", "59Peg",
    "mu.For", "HR1075", "HR2160", "alfPic", "HD71043", "KCen", "HR4977", "tauCen",
    "c02Cen", "HD129422", "HR5494", "HD131625", "lamNor", "HD154310", "alfCrA", "HR8366",
    "mu.PsA"
]
TELLURIC_TARGETS['NIRPS_HE'] = [
    "HR9098", "HR806", "HR1903", "HR3117", "HR3314", "HR4023", "HR4467", "HR4468",
    "HR4889", "HR5671", "HR6743", "HR7590", "HR8709", "HR875", "HR3131", "HR5107",
    "HR7830", "15_PEG", "ZETVIR"
]
# NIRPS-HE and NIRPS-HA have the same telluric targets
TELLURIC_TARGETS['NIRPS_HA'] = list(TELLURIC_TARGETS['NIRPS_HE'])
# -----------------------------------------------------------------------------
# set the list of observation directories to check for each instrument
OBSERVATION_DIRS = dict()
OBSERVATION_DIRS['SPIROU'] = ["2020-05-14", "2020-06-30", "2020-07-29",
                              "2020-07-30", "2020-08-01", "2020-08-02",
                              "2020-08-03", "2020-08-04", "2020-08-10",
                              "2020-09-23", "2020-10-07", "2020-11-01"]
OBSERVATION_DIRS['NIRPS_HE'] = None
OBSERVATION_DIRS['NIRPS_HA'] = None
# -----------------------------------------------------------------------------
# Set the reference observation directory for each instrument (this night must
#   be included)
REF_OBSERVATION_DIRS = dict()
REF_OBSERVATION_DIRS['SPIROU'] = "2020-08-31"
REF_OBSERVATION_DIRS['NIRPS_HE'] = "2022-11-24"
REF_OBSERVATION_DIRS['NIRPS_HA'] = "2022-11-24"
# -----------------------------------------------------------------------------
# Set the maximum number of observation directories we want for each instrument
MAX_OBSERVATION_DIRS = dict()
MAX_OBSERVATION_DIRS['SPIROU'] = 5
MAX_OBSERVATION_DIRS['NIRPS_HE'] = 5
MAX_OBSERVATION_DIRS['NIRPS_HA'] = 5

# Set the maximum number of science targets we want for each instrument
MAX_SCI_OBJS = dict()
MAX_SCI_OBJS['SPIROU'] = 1
MAX_SCI_OBJS['NIRPS_HE'] = 1
MAX_SCI_OBJS['NIRPS_HA'] = 1

# Set the maximum number of telluric stars we want for each instrument
MAX_TELLU_OBJS = dict()
MAX_TELLU_OBJS['SPIROU'] = 3
MAX_TELLU_OBJS['NIRPS_HE'] = 3
MAX_TELLU_OBJS['NIRPS_HA'] = 3
# -----------------------------------------------------------------------------

REQCALS = dict(SPIROU=[], NIRPS_HE=[], NIRPS_HA=[])
REQCALS['SPIROU'].append(('DARK_DARK_TEL', 2))
REQCALS['SPIROU'].append(('DARK_DARK_INT', 2))
REQCALS['SPIROU'].append(('DARK_FLAT', 2))
REQCALS['SPIROU'].append(('DARK_FP', 2))
REQCALS['SPIROU'].append(('FLAT_FLAT', 2))
REQCALS['SPIROU'].append(('FLAT_DARK', 2))
REQCALS['SPIROU'].append(('FP_FP', 2))
REQCALS['SPIROU'].append(('HCONE_HCONE', 2))

REQCALS['NIRPS_HE'].append(('DARK_DARK', 1))
REQCALS['NIRPS_HE'].append(('DARK_FLAT', 1))
REQCALS['NIRPS_HE'].append(('DARK_FP', 1))
REQCALS['NIRPS_HE'].append(('FLAT_FLAT', 1))
REQCALS['NIRPS_HE'].append(('FLAT_DARK', 1))
REQCALS['NIRPS_HE'].append(('FP_FP', 1))
REQCALS['NIRPS_HE'].append(('HCONE_HCONE', 1))

REQCALS['NIPRS_HA'] = list(REQCALS['NIRPS_HE'])

# Do not ask for directories (only for use after a test)
SELECTED_OBS_DIRS = dict()
SELECTED_OBS_DIRS['SPIROU'] = ['2020-08-31', '2020-08-10', '2020-07-30',
                               '2020-05-14', '2020-10-07']
SELECTED_OBS_DIRS['NIRPS_HE'] = []
SELECTED_OBS_DIRS['NIRPS_HA'] = []


# =============================================================================
# Define functions
# =============================================================================
def obs_from_raw_dir(params) -> List[str]:
    # get the raw directory
    raw_dir = params['DRS_DATA_RAW']
    # get everything in the raw directory
    raw_obs_dirs = os.listdir(raw_dir)
    # filter these to make sure they are directories
    obs_dirs = []
    for raw_obs_dir in raw_obs_dirs:
        raw_obs_path = os.path.join(raw_dir, raw_obs_dir)
        if os.path.isdir(raw_obs_path):
            obs_dirs.append(raw_obs_dir)
    # return obs_dirs
    return obs_dirs


def clean_objs(params, pconst, database, raw_objs) -> List[str]:
    # print progress
    WLOG(params, '', 'Cleaning {0} objects'.format(len(raw_objs)))
    # clean object names to match APEROs object names
    allowed_objs = []
    for raw_obj in tqdm(raw_objs):
        correct_objname, found = database.find_objname(pconst, raw_obj)
        if found > 0:
            allowed_objs.append(correct_objname)
    return allowed_objs


def get_valid_obs_dirs(params, obs_dirs: List[str]) -> List[str]:
    # get the instrument
    instrument = params['INSTRUMENT']
    # deal with None
    if obs_dirs is None:
        obs_dirs = obs_from_raw_dir(params)
    # get the database
    findexdb = drs_database.FileIndexDatabase(params)
    # load the database
    findexdb.load_db()
    # set up condition (raw data only)
    condition = 'BLOCK_KIND="raw"'
    # add the observations
    subconditions = []
    for observation in obs_dirs:
        subconditions.append(f'OBS_DIR="{observation}"')
    # add the reference observation to the subconditions
    subconditions.append(f'OBS_DIR="{REF_OBSERVATION_DIRS[instrument]}"')
    # add to condition
    condition += ' AND ({0})'.format(' OR '.join(subconditions))
    # get a table for this condition
    table = findexdb.get_entries('ABSPATH, OBS_DIR, KW_DPRTYPE, KW_OBJNAME',
                                 condition=condition)
    # list of unique observation directories (we may not have all those given)
    uobsdirs = list(set(table['OBS_DIR']))
    # storage for valid observation directories
    valid_obs_dirs = []
    # loop around unique observation directories
    for uobsdir in uobsdirs:
        # get the unique table for this obsdir
        utable = table[table['OBS_DIR'] == uobsdir]
        # get the unique calibrations
        cals = list(utable['KW_DPRTYPE'])
        # assume directory has passed
        passed = True
        # loop around required calibrations and check if they are present
        #   if not our directory has not passed
        for cal in REQCALS[instrument]:
            calname, calnum = cal
            if calname not in cals:
                passed = False
                break
            if cals.count(calname) < calnum:
                passed = False
                break
        # if passed all requirements then pass it
        if passed:
            valid_obs_dirs.append(uobsdir)

    # double check that the ref observation directory is still in the list
    if REF_OBSERVATION_DIRS[instrument] not in valid_obs_dirs:
        emsg = 'Reference observation directory = {0} not valid'
        eargs = [REF_OBSERVATION_DIRS[instrument]]
        WLOG(params, 'error', emsg.format(*eargs))

    # return a list of valid observation directories
    return valid_obs_dirs


def display_obs_dirs(params, obs_dirs, count_sci_obj, count_tellu_obj,
                     best_sci_objs, best_tellu_objs) -> List[str]:
    # print progress
    WLOG(params, '', 'Listing number of files for each object')
    # storage for options (these are the numbers the user is choosing from)
    options = []
    # need a counter to tell user which option they are selecting
    option_counter = 0
    # loop around observation directories
    for o_it, obsdir in enumerate(obs_dirs):
        # storage for targets found
        printouts = []
        # storage for number of science and telluric objects found in each
        #   observation directory
        sci_obj_found = 0
        tellu_obj_found = 0
        # loop around science directories
        for sci_obj in best_sci_objs:
            if count_sci_obj[obsdir][sci_obj] == 0:
                continue
            sci_obj_found += count_sci_obj[obsdir][sci_obj]
            if sci_obj_found > 0:
                wmsg1 = '\t{0} = {1}'
                wargs1 = [sci_obj, count_sci_obj[obsdir][sci_obj]]
                printouts.append(wmsg1.format(*wargs1))
        # loop around telluric directories
        for tellu_obj in best_tellu_objs:
            if count_tellu_obj[obsdir][tellu_obj] == 0:
                continue
            tellu_obj_found += count_tellu_obj[obsdir][tellu_obj]
            if tellu_obj_found > 0:
                wmsg2 = '\t{0} = {1}'
                wargs2 = [tellu_obj, count_tellu_obj[obsdir][tellu_obj]]
                printouts.append(wmsg2.format(*wargs2))
        # do not show observation directory if no science or telluric objects
        #   present
        if sci_obj_found == 0 or tellu_obj_found == 0:
            continue
        else:
            # print progress
            wmsg = '{0}: OBSDIR = {1}'
            wargs = [option_counter + 1, obsdir]
            WLOG(params, 'info', wmsg.format(*wargs))

            for printout in printouts:
                WLOG(params, '', printout)

            options.append(obsdir)
            option_counter += 1

    return options


def count_files(params, obs_dirs, condition, sci_objs, tellu_objs):
    # print progress
    wmsg = 'Counting files for {0} observation directories'
    wargs = [len(obs_dirs)]
    WLOG(params, '', wmsg.format(*wargs))
    # get the database
    findexdb = drs_database.FileIndexDatabase(params)
    # load the database
    findexdb.load_db()
    # storage for counts
    count_total_sci = dict()
    count_total_tellu = dict()
    count_sci_obj = dict()
    count_tellu_obj = dict()
    # loop around each observation directory
    for obsdir in tqdm(obs_dirs):
        # set up condition to add this observation directory
        condition1 = condition + f' AND OBS_DIR="{obsdir}"'
        # add obsdir to count dict
        count_sci_obj[obsdir] = dict()
        count_tellu_obj[obsdir] = dict()
        # count science objects
        for sci_obj in sci_objs:
            condition2 = condition1 + f' AND KW_OBJNAME="{sci_obj}"'
            count = findexdb.database.count('*', condition=condition2)
            count_sci_obj[obsdir][sci_obj] = count
            if count == 0:
                continue
            if sci_obj in count_total_sci:
                count_total_sci[sci_obj] += count
            else:
                count_total_sci[sci_obj] = count
        # count telluric objects
        for tellu_obj in tellu_objs:
            condition2 = condition1 + f' AND KW_OBJNAME="{tellu_obj}"'
            count = findexdb.database.count('*', condition=condition2)
            count_tellu_obj[obsdir][tellu_obj] = count
            if count == 0:
                continue
            if tellu_obj in count_total_tellu:
                count_total_tellu[tellu_obj] += 1
            else:
                count_total_tellu[tellu_obj] = 1

    return count_total_sci, count_total_tellu, count_sci_obj, count_tellu_obj


def select_obs_dirs(params, obs_dirs, count_sci_obj, count_tellu_obj,
                    best_sci_objs, best_tellu_objs):
    # get instrument
    instrument = params['INSTRUMENT']
    # sort obs_dirs alphabetically
    out_obs_dirs = np.sort(obs_dirs)
    # copy this to remove from list
    left_obs_dirs = list(out_obs_dirs)
    # storage of selected observation directories
    selected = []

    while len(selected) < MAX_OBSERVATION_DIRS[instrument]:

        # print message
        msg = 'Please select the best observation directory ({0} remaining)'
        margs = [MAX_OBSERVATION_DIRS[instrument] - len(selected)]
        WLOG(params, '', msg.format(*margs), colour='magenta')

        options = display_obs_dirs(params, left_obs_dirs, count_sci_obj,
                                   count_tellu_obj, best_sci_objs,
                                   best_tellu_objs)

        qmsg = ('Enter the number of the observation directory you want to '
                'select [1 - {0}]: ')
        qargs = [len(options)]
        uinput = input(qmsg.format(*qargs))

        try:
            uinput = int(uinput) - 1
            selected.append(str(options[uinput]))
            left_obs_dirs.remove(options[uinput])
        except Exception as _:
            wmsg = 'Input = {0} is not a valid integer between 1 and {1}'
            wargs = [uinput, len(left_obs_dirs)]
            WLOG(params, 'warning', wmsg.format(*wargs))
            continue

    # deal with reference observation directory not selected (add it in)
    if REF_OBSERVATION_DIRS[instrument] not in selected:
        selected.append(REF_OBSERVATION_DIRS[instrument])

    return selected


class RawFile:
    def __init__(self, kind, abspath, obsdir, filename, objname, dprtype,
                 mjdmid):
        self.kind = kind
        self.abspath = abspath
        self.obsdir = obsdir
        self.filename = filename
        self.objname = objname
        self.dprtype = dprtype
        self.mjdmid = mjdmid

    def __str__(self):
        return 'RawFile[{0}]={1}'.format(self.dprtype, self.filename)

    def __repr__(self):
        return self.__str__()


class ObsDir:
    def __init__(self, obsdir):
        self.obsdir = obsdir
        self.files = []
        self.sci_objs = []
        self.tellu_objs = []
        self.calibs = []

    def add_file(self, rawfile):
        self.files.append(rawfile)
        if rawfile.kind == 'sci':
            self.sci_objs.append(rawfile)
        if rawfile.kind == 'tellu':
            self.tellu_objs.append(rawfile)
        if rawfile.kind == 'calib':
            self.calibs.append(rawfile)

    def __str__(self):
        return f'ObsDir[{self.obsdir}]'

    def __repr__(self):
        return self.obsdir

    def cut_cals(self, reqcals):

        all_calibs = []

        for reqcal in reqcals:

            udprtype, number_of_cals = reqcal

            calibs = list(filter(lambda x: x.dprtype == udprtype, self.calibs))

            # get a sort mask based on mjdmid (time)
            argsort = np.array(list(map(lambda x: x.mjdmid, calibs))).argsort()
            # sort by mjdmid (time)
            prev_calibs = np.array(calibs)[argsort]
            # set up a new set of calibrations
            keep_calibs = []
            # go from each end (first then last) and add to calibs
            while len(keep_calibs) < number_of_cals and len(prev_calibs) > 0:
                if len(keep_calibs) % 2 == 0:
                    keep_calibs.append(prev_calibs[0])
                    prev_calibs = prev_calibs[1:]
                else:
                    keep_calibs.append(prev_calibs[-1])
                    prev_calibs = prev_calibs[:-1]
            # push into full list
            all_calibs += keep_calibs
        # push back into calibs
        self.calibs = all_calibs

    def get_copy_paths(self, outpath):

        inpaths = []
        outpaths = []

        for sci_obj in self.sci_objs:
            inpaths.append(sci_obj.abspath)
            outpaths.append(os.path.join(outpath, self.obsdir, sci_obj.filename))
        for tellu_obj in self.tellu_objs:
            inpaths.append(tellu_obj.abspath)
            outpaths.append(os.path.join(outpath, self.obsdir, tellu_obj.filename))
        for calib in self.calibs:
            inpaths.append(calib.abspath)
            outpaths.append(os.path.join(outpath, self.obsdir, calib.filename))

        return inpaths, outpaths


def get_selected_files(params, selected_obsdirs, best_sci_objs, best_tellu_objs,
                       reqcals) -> List[ObsDir]:
    # get the database
    findexdb = drs_database.FileIndexDatabase(params)
    # load the database
    findexdb.load_db()
    # storage for file list
    file_list = []
    # loop around selected observation directories
    for obsdir in selected_obsdirs:
        # create class
        obsdir_inst = ObsDir(obsdir)
        # set up condition to add this observation directory
        condition = 'BLOCK_KIND="raw" AND OBS_DIR="{0}"'.format(obsdir)
        # get the table for this observation directory
        cols = ('ABSPATH, OBS_DIR, FILENAME, KW_OBJNAME, KW_DPRTYPE, '
                'KW_MID_OBS_TIME')
        table = findexdb.get_entries(cols, condition=condition)
        # loop around the table and add files to file_list
        for row in range(len(table)):
            abspath, obsdir, filename, objname, dprtype, mjdmid = table.iloc[row]
            # assume we don't want file
            kind = None
            # check if this is a science object
            if objname in best_sci_objs:
                kind = 'sci'
            # check if this is a telluric object
            if objname in best_tellu_objs:
                kind = 'tellu'
            # check if this is a required calibration
            for reqcal in reqcals:
                if dprtype == reqcal[0]:
                    kind = 'calib'
                    break
            # add file if we have a kind (i.e. we want the file)
            if kind is not None:
                rawfile = RawFile(kind, abspath, obsdir, filename, objname,
                                  dprtype, mjdmid)
                obsdir_inst.add_file(rawfile)
        # cut down the number of calibrations files to the required number of
        #   calibraitons
        obsdir_inst.cut_cals(reqcals)
        # push obsdir class into file_list
        file_list.append(obsdir_inst)
    return file_list


def main(params: ParamDict):
    # load pconst
    pconst = constants.pload()
    # get the instrument
    instrument = params['INSTRUMENT']

    # get the object database
    objdb = drs_database.AstrometricDatabase(params)
    # load the object database
    objdb.load_db()

    # ----------------------------------------------------------------------
    # step 1: get a valid list of observation directories
    # ----------------------------------------------------------------------
    # rule out any observation directories that do not have required calibrations
    obs_dirs = get_valid_obs_dirs(params, OBSERVATION_DIRS[instrument])
    # set up condition (raw data only)
    condition = 'BLOCK_KIND="raw"'

    # -------------------------------------------------------------------------
    # get a list of allowed object names
    raw_sci = SCIENCE_TARGETS[instrument]
    raw_tellu = TELLURIC_TARGETS[instrument]
    # clean object names (to match APERO)
    sci_objs = clean_objs(params, pconst, objdb, raw_sci)
    tellu_objs = clean_objs(params, pconst, objdb, raw_tellu)

    # -------------------------------------------------------------------------
    # We need statistics on how many files we have for each directory
    # -------------------------------------------------------------------------
    cout = count_files(params, obs_dirs, condition, sci_objs, tellu_objs)
    count_total_sci, count_total_tellu, count_sci_obj, count_tellu_obj = cout

    # -------------------------------------------------------------------------
    # find the best science target (highest number of files in count_total)
    # -------------------------------------------------------------------------
    # sort the count_sci_obj by the total number of files
    sci_values = list(count_total_sci.values())
    sci_keys = list(count_total_sci.keys())
    sci_sort = np.argsort(sci_values)[-MAX_SCI_OBJS[instrument]:]
    best_sci_objs = np.array(sci_keys)[sci_sort]

    # -------------------------------------------------------------------------
    # find the 5 best telluric targets
    # -------------------------------------------------------------------------
    # sort the count_tellu_obj by the total number of files
    tellu_values = list(count_total_tellu.values())
    tellu_keys = list(count_total_tellu.keys())
    tellu_sort = np.argsort(tellu_values)[-MAX_TELLU_OBJS[instrument]:]
    best_tellu_objs = np.array(tellu_keys)[tellu_sort]

    # -------------------------------------------------------------------------
    # Ask user to select best observation directories
    # -------------------------------------------------------------------------
    if len(SELECTED_OBS_DIRS[instrument]) == 0:

        selected_obsdirs = select_obs_dirs(params, obs_dirs, count_sci_obj,
                                           count_tellu_obj, best_sci_objs,
                                           best_tellu_objs)
    else:
        selected_obsdirs = SELECTED_OBS_DIRS[instrument]

    # -------------------------------------------------------------------------
    # get files for this night
    #   1. best science targets
    #   2. best telluric targets
    #   3. calibrations (cut down using REQCAL)
    #   for each observation directory in selected
    # -------------------------------------------------------------------------
    file_list = get_selected_files(params, selected_obsdirs, best_sci_objs,
                                   best_tellu_objs, REQCALS[instrument])

    # print which nights were selected
    WLOG(params, '', 'Selected observation directories:')
    for obsdir_inst in file_list:
        WLOG(params, '', '\t - {0}'.format(obsdir_inst.obsdir))

    # -------------------------------------------------------------------------
    # copy files to destination
    # -------------------------------------------------------------------------
    outdir = OUTPUTDIR[instrument]
    # loop around each observation directory
    for obsdir_inst in file_list:
        # make directory if it doesn't exist
        if not TEST:
            obsdir_path = os.path.join(outdir, obsdir_inst.obsdir)
            if not os.path.exists(obsdir_path):
                msg = 'Creating directory {0}'
                margs = [obsdir_path]
                WLOG(params, '', msg.format(*margs))
                os.makedirs(obsdir_path)
        # construct inpaths and outpaths
        inpaths, outpaths = obsdir_inst.get_copy_paths(outdir)
        # copy files
        for inpath, outpath in zip(inpaths, outpaths):
            msg = 'Copying {0} to {1}'
            margs = [inpath, outpath]
            WLOG(params, '', msg.format(*margs))
            # copy file
            if not TEST:
                shutil.copy(inpath, outpath)

    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load apero parameters
    _params = constants.load()
    _params.set('PID', 'Unknown')
    # run the main function
    main(_params)


# =============================================================================
# End of code
# =============================================================================
