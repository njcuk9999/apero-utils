#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-06-11

@author: cook
"""
import os
import shutil
from apero.core import constants
from apero.core.core import drs_database


# =============================================================================
# Define variables
# =============================================================================
# list all object names that should be found
kw_objnames = ['HD189733', 'WASP-127']
# kw_objnames = ['TOI 1449']
# list all DPRTYPEs that should be found
kw_dprtypes = ['OBJ_FP', 'OBJ_DARK', 'POLAR_FP', 'POLAR_DARK']
# list all drs output types that should be found
kw_outputs = ['EXT_E2DS_FF']
kw_outputs = ['EXT_E2DS_FF', 'TELLU_OBJ', 'SC1D_V_FILE', 'SC1D_W_FILE']
kw_outputs = ['TELLU_OBJ', 'TELLU_RECON']
# list all fibers required
kw_fibers = ['AB']
# outpath (where we should copy these to)
OUTPATH = '/home/boucher/spirou/planetModels/HD_189733_b/full_210523'
#OUTPATH = '/spirou2/full_210523_objs'
# copy (debug set to False)
COPY = True

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # load the constants
    print('Starting APERO...')
    params = constants.load()
    pconst = constants.pload()
    # load index database
    print('Loading database...')
    indexdb = drs_database.IndexDatabase(params)
    indexdb.load_db()

    # -------------------------------------------------------------------------
    # create master condition
    master_condition = ''
    # ---------------------------------------------------------------------
    # dprtypes
    # ---------------------------------------------------------------------
    if kw_dprtypes is not None:
        subconditions = []
        # loop around object names
        for kw_dprtype in kw_dprtypes:
            subcondition = '(KW_DPRTYPE="{0}")'.format(kw_dprtype)
            subconditions.append(subcondition)
        # add to condition
        master_condition += '({0})'.format(' OR '.join(subconditions))
    # ---------------------------------------------------------------------
    # outputs
    # ---------------------------------------------------------------------
    if kw_outputs is not None:
        subconditions = []
        # loop around object names
        for kw_output in kw_outputs:
            subcondition = '(KW_OUTPUT="{0}")'.format(kw_output)
            subconditions.append(subcondition)
        # add to condition
        master_condition += ' AND ({0})'.format(' OR '.join(subconditions))
    # ---------------------------------------------------------------------
    # fibers
    # ---------------------------------------------------------------------
    if kw_fibers is not None:
        subconditions = []
        # loop around object names
        for kw_fiber in kw_fibers:
            subcondition = '(KW_FIBER="{0}")'.format(kw_fiber)
            subconditions.append(subcondition)
        # add to condition
        master_condition += ' AND ({0})'.format(' OR '.join(subconditions))
    # -------------------------------------------------------------------------
    # separate list for each object name
    # -------------------------------------------------------------------------
    # storage of inpaths
    database_inpaths = dict()
    # loop around input object names
    for kw_objname in kw_objnames:
        # clean object name (as best we can)
        clean_obj_name = pconst.DRS_OBJ_NAME(kw_objname)
        print('Processing KW_OBJNAME={0}'.format(clean_obj_name))
        # write condition for this object
        obj_condition = '(KW_OBJNAME="{0}")'.format(clean_obj_name)
        # deal with having a master condition
        if len(master_condition) > 0:
            condition = '{0} AND {1}'.format(obj_condition, master_condition)
        else:
            condition = str(obj_condition)
        # get inpaths
        inpaths = indexdb.get_entries('ABSPATH', condition=condition)
        # load into file storage
        if len(inpaths) > 0:
            print('\tFound {0} entries'.format(len(inpaths)))
            database_inpaths[clean_obj_name] = inpaths
        else:
            print('\tFound no entries')
    # -------------------------------------------------------------------------
    # Now get outpaths (if infile exists)
    # -------------------------------------------------------------------------
    # storage of inpaths/outpaths
    all_inpaths = dict()
    all_outpaths = dict()
    # loop around objects with files
    for objname in database_inpaths:
        # output directory for objname
        outdir = os.path.join(OUTPATH, objname)
        # print progress
        print('Adding outpaths for KW_OBJNAME={0}'.format(objname))
        # add object name to storage
        all_inpaths[objname] = []
        all_outpaths[objname] = []
        # loop around all files for this object
        for filename in database_inpaths[objname]:
            # if object exists
            if os.path.exists(filename):
                # get paths
                inpath = filename
                basename = os.path.basename(filename)
                outpath = os.path.join(outdir, basename)
                # add to storage
                all_inpaths[objname].append(inpath)
                all_outpaths[objname].append(outpath)
        # make a directory for this object (if it doesn't exist)
        if len(all_outpaths[objname]) != 0:
            # print progress
            print('\tAdded {0} outpaths'.format(len(all_outpaths[objname])))
            # create output directory if it doesn't exist
            if not os.path.exists(outdir) and COPY:
                os.mkdir(outdir)
    # -------------------------------------------------------------------------
    # Copy files
    # -------------------------------------------------------------------------
    for objname in all_inpaths:
        print()
        print('=' * 50)
        print('COPY OBJNAME={0}'.format(objname))
        print('=' * 50)
        print()
        # loop around files
        for row in range(len(all_inpaths[objname])):
            # get in and out path
            inpath = all_inpaths[objname][row]
            outpath = all_outpaths[objname][row]
            # print string
            copyargs = [row + 1, len(all_inpaths[objname]), outpath]
            copystr = '[{0}/{1}] --> {2}'.format(*copyargs)
            # print copy string
            print(copystr)
            # copy
            if COPY:
                shutil.copy(inpath, outpath)

# =============================================================================
# End of code
# =============================================================================
