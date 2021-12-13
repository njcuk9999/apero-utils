#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-11-10 at 16:04

@author: cook
"""
from astroquery.utils.tap.core import TapPlus
from astroquery.simbad import Simbad
import numpy as np
from typing import Dict, List, Tuple
import warnings

from apero import lang
from apero.base import base
from apero.core import constants
from apero.core.core import drs_log
from apero.core.core import drs_database
from apero.tools.module.database import manage_databases

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'apero_astrometrics.py'
__INSTRUMENT__ = 'None'
__PACKAGE__ = base.__PACKAGE__
__version__ = base.__version__
__author__ = base.__author__
__date__ = base.__date__
__release__ = base.__release__
# get text entry instance
textentry = lang.textentry
# Get Logging function
WLOG = drs_log.wlog
# get the databases
IndexDatabase = drs_database.IndexDatabase
ObjectDatabase = drs_database.ObjectDatabase
# simbad additional columns
SIMBAD_COLUMNS = ['ids']


# =============================================================================
# Define functions
# =============================================================================
def query_object(rawobjname):
    # set up the TapPlus class
    # simbad = TapPlus(url=URL)

    # execute the job
    # sjob = simbad.launch_job(QUERY.format(rawobjname))

    # get the results
    # table = sjob.get_results()

    # get results
    with warnings.catch_warnings(record=True) as _:
        # add ids column
        for simbad_column in SIMBAD_COLUMNS:
            Simbad.add_votable_fields(simbad_column)
        # query simbad
        table = Simbad.query_object(rawobjname)

    return 0



def query_database(params, rawobjnames: List[str] ) -> List[str]:
    """
    Find all objects in the object database and return list of unfound
    objects

    :param params:
    :param rawobjnames:
    :return:
    """
    # ---------------------------------------------------------------------
    # get psuedo constants
    pconst = constants.pload()
    # ---------------------------------------------------------------------
    # Update the object database (recommended only for full reprocessing)
    # check that we have entries in the object database
    manage_databases.object_db_populated(params)
    # update the database if required
    if params['UPDATE_OBJ_DATABASE']:
        # log progress
        WLOG(params, '', textentry('40-503-00039'))
        # update database
        manage_databases.update_object_database(params, log=False)
    # ---------------------------------------------------------------------
    # print progress
    WLOG(params, '', 'Searching local object database for object names...')
    # load the object database after updating
    objdbm = ObjectDatabase(params)
    objdbm.load_db()
    # storage for output - assume none are found
    unfound = []
    # loop around objects and find them in database
    for rawobjname in rawobjnames:
        # clean the object
        objname = pconst.DRS_OBJ_NAME(rawobjname)
        # find correct name in the database (via objname or aliases)
        correct_objname, found = objdbm.find_objname(pconst, objname)
        # deal with found / not found
        if found:
            msg = '\tObject: "{0}" found in database as "{1}"'
            margs = [rawobjname, correct_objname]
            WLOG(params, '', msg.format(*margs))
        else:
            msg = '\tObject: "{0}" not found in database'
            margs = [rawobjname]
            WLOG(params, '', msg.format(*margs), colour='yellow')
            # add to unfound list
            unfound.append(rawobjname)
    # return the entry names and the found dictionary
    return unfound


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get params
    params = constants.load()
    params.set('PID', 'None')
    params.set('UPDATE_OBJ_DATABASE', False)

    rawobjs = ['Gl699', 'Trappist1', 'Neil']

    # ----------------------------------------------------------------------
    # step 1: Is object in database?
    # ----------------------------------------------------------------------
    # query local object database
    unfound_objs = query_database(params, rawobjs)

    # stop here if all objects found



    print('stop')
    # ----------------------------------------------------------------------
    # step 2: if not in database see if we can find it in simbad
    # _ = query_object('Gl699')


    # ----------------------------------------------------------------------
    # step 3: add to pending database if submitted

# =============================================================================
# End of code
# =============================================================================
