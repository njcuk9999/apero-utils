#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-02 at 11:47

@author: cook
"""

import sys
from typing import List

from apero import lang
from apero.base import base
from apero.core import constants
from apero.core.core import drs_log
from apero.core.core import drs_database
from apero.core.instruments.default import pseudo_const

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
# Define the required science targets for a mini dataset
SCIENCE_TARGETS = dict()
SCIENCE_TARGETS['SPIROU'] = ['GL699']
SCIENCE_TARGETS['NIRPS_HE'] = ['PROXIMA']
SCIENCE_TARGETS['NIRPS_HA'] = ['PROXIMA']
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
OBSERVATION_DIRS['SPIROU'] = [
    "2020-05-14", "2020-06-30", "2020-07-29", "2020-07-30", "2020-08-01",
    "2020-08-02", "2020-08-03", "2020-08-04", "2020-08-10",
    "2020-09-23", "2020-10-07", "2020-11-01"
]
OBSERVATION_DIRS['NIRPS_HE'] = []
OBSERVATION_DIRS['NIRPS_HA'] = []
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
# -----------------------------------------------------------------------------
# set the required number of calibrations for each instrument
REQUIRED_CALIBRATIONS = dict()
REQUIRED_CALIBRATIONS['SPIROU'] = 2
REQUIRED_CALIBRATIONS['NIRPS_HE'] = 1
REQUIRED_CALIBRATIONS['NIRPS_HA'] = 1


# =============================================================================
# Define functions
# =============================================================================
def clean_objs(pconst, database, raw_objs) -> List[str]:
    # clean object names to match APEROs object names
    allowed_objs = []
    for raw_obj in raw_objs:
        correct_objname, found = database.find_objname(pconst, raw_obj)
        if found > 0:
            allowed_objs.append(correct_objname)
    return allowed_objs


def main(params: ParamDict):
    # load pconst
    pconst = constants.pload()
    # get the instrument
    instrument = params['INSTRUMENT']
    # get the database
    findexdb = drs_database.FileIndexDatabase(params)
    # load the database
    findexdb.load_db()

    # get the object database
    objdb = drs_database.AstrometricDatabase(params)
    # load the object database
    objdb.load_db()

    # ----------------------------------------------------------------------
    # step 1: get a table of all raw files that match the instrument
    # ----------------------------------------------------------------------
    # set up condition (raw data only)
    condition = 'BLOCK_KIND="raw"'
    # add the observations
    subconditions = []
    for observation in OBSERVATION_DIRS[instrument]:
        subconditions.append(f'OBS_DIR="{observation}"')
    # add the reference observation to the subconditions
    subconditions.append(f'OBS_DIR="{REF_OBSERVATION_DIRS[instrument]}"')
    # add to condition
    condition += ' AND ({0})'.format(' OR '.join(subconditions))


    # get a table for this condition
    # table = findexdb.get_entries('ABSPATH, OBS_DIR, KW_DPRTYPE, KW_OBJNAME',
    #                              condition=condition)

    # -------------------------------------------------------------------------
    # get a list of allowed object names
    raw_sci = SCIENCE_TARGETS[instrument]
    raw_tellu = TELLURIC_TARGETS[instrument]
    # clean object names (to match APERO)
    sci_objs = clean_objs(pconst, objdb, raw_sci)
    tellu_objs = clean_objs(pconst, objdb, raw_tellu)
    # -------------------------------------------------------------------------
    # We need statistics on how many files we have for each directory
    # -------------------------------------------------------------------------
    for obsdir in OBSERVATION_DIRS[instrument]:

        condition1 = condition + f' AND OBS_DIR="{obsdir}"'

        count_total = dict()
        count_sci_obj = dict()
        count_tellu_obj = dict()
        # count science objects
        for sci_obj in sci_objs:
            condition2 = condition1 + f' AND KW_OBJNAME="{sci_obj}"'
            count = findexdb.database.count('*', condition=condition2)
            count_sci_obj[sci_obj] = count

            if sci_obj in count_total:
                count_total[sci_obj] += count
            else:
                count_total[sci_obj] = count
        # count telluric objects
        for tellu_obj in tellu_objs:
            condition2 = condition1 + f' AND KW_OBJNAME="{tellu_obj}"'
            count = findexdb.database.count('*', condition=condition2)
            count_tellu_obj[tellu_obj] = count

            if tellu_obj in count_total:
                count_total[tellu_obj] += count
            else:
                count_total[tellu_obj] = count


    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load apero parameters
    _params = constants.load()
    # run the main function
    main(_params)


# =============================================================================
# End of code
# =============================================================================
