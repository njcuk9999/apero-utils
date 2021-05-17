#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-05-10

@author: cook
"""
from astropy.table import Table
from astropy.io import fits
from pathlib import Path

from apero.core import constants
from apero.core.core import drs_file


# =============================================================================
# Define variables
# =============================================================================
filename = '/data/spirou/data/07000/calibDB/8299B45D9D_pp_flat_AB.fits'
filename = '/data/spirou/data/07000/tmp/2019-04-23/2401201f_pp.fits'

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # load params
    params = constants.load()
    # load pseduo constants
    pconst = constants.pload()
    # -------------------------------------------------------------------------
    # load param table
    ptable = Table.read(filename, hdu='PARAM_TABLE')
    # load the header
    header = fits.getheader(filename)
    # get drs path
    drspath = drs_file.DrsPath(params, filename)
    # get path instance
    filepath = Path(filename)
    # get block for this path
    block = drspath.get_block()
    lastmodified = filepath.stat().st_mtime

    # -------------------------------------------------------------------------
    # make log entry
    # -------------------------------------------------------------------------
    # log entry mask
    logmask = ptable['KIND'] == 'rlog'
    # push into a dictionary (for easy access)
    logdict = dict()
    for row, key in enumerate(ptable[logmask]['NAME']):
        logdict[key] = ptable[logmask]['VALUE'][row]
    # get log keys and types
    logcols, _ = pconst.LOG_DB_COLUMNS()
    # loop around log keys and add them to values
    logvalues = []
    for logkey in logcols:
        # construct keys
        key = 'rlog.{0}'.format(logkey.lower())
        # get value
        logvalue = logdict.get(key, 'NULL')
        # append value to values
        logvalues.append(logvalue)
    # TODO: now just need to add this to database

    # -------------------------------------------------------------------------
    # make index entry
    # -------------------------------------------------------------------------
    # only index if we are supposed to index this block kind
    if block.indexing:
        # get index columns
        indexcols, _, _ = pconst.INDEX_DB_COLUMNS()
        indexheaderkeys, _ = pconst.INDEX_HEADER_KEYS()
        # store index values
        indexvalues = []
        # abspath
        indexvalues.append(drspath.abspath)
        # get obs dir
        indexvalues.append(drspath.obs_dir)
        # add block kind
        indexvalues.append(drspath.block_kind)
        # get filename
        indexvalues.append(drspath.basename)
        # get last modified time
        indexvalues.append(lastmodified)
        # get recipe name
        indexvalues.append(logdict['rlog.recipe'])
        # get runstring
        indexvalues.append(logdict['rlog.runstring'])
        # loop around index header keys
        for indexcol in indexheaderkeys:
            # get key from params
            if indexcol in params:
                drskey = params[indexcol][0]
            else:
                drskey = indexcol

            # get these keys from header
            hvalue = header.get(drskey, 'NULL')
            # add header key to index values
            indexvalues.append(hvalue)
        # add used and rawfix
        indexvalues.append(1)
        indexvalues.append(1)

    # TODO: now just need to add this to database

    # -------------------------------------------------------------------------
    # make calib entry
    # -------------------------------------------------------------------------
    if block.name == 'calib':
        pass

    # -------------------------------------------------------------------------
    # make tellu entry
    # -------------------------------------------------------------------------
    if block.name == 'tellu':
        pass

# =============================================================================
# End of code
# =============================================================================
